#!/usr/bin/env python3
"""
Diagnostic script to test AWS detection methods.

Run this on your JupyterHub (or any environment) to see which
detection methods succeed. Copy-paste into a notebook cell or
run from a terminal:

    python test_aws_detection_diagnostic.py
"""

import os
import sys


def test_imdsv2():
    """Test 1: EC2 Instance Metadata Service v2 (current method)."""
    print("=" * 60)
    print("TEST 1: EC2 IMDSv2 metadata service")
    print("=" * 60)
    try:
        import requests
        # IMDSv2 token request
        token_resp = requests.put(
            'http://169.254.169.254/latest/api/token',
            headers={'X-aws-ec2-metadata-token-ttl-seconds': '60'},
            timeout=2)
        token_resp.raise_for_status()
        token = token_resp.text
        print(f"  Token obtained: {token[:20]}...")

        meta_resp = requests.get(
            'http://169.254.169.254/latest/meta-data/instance-id',
            headers={'X-aws-ec2-metadata-token': token},
            timeout=2)
        meta_resp.raise_for_status()
        print(f"  Instance ID: {meta_resp.text}")
        print("  RESULT: DETECTED as AWS EC2")
        return True
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        print("  RESULT: NOT detected via IMDSv2")
        return False


def test_imdsv1():
    """Test 2: EC2 Instance Metadata Service v1 (legacy, no token)."""
    print()
    print("=" * 60)
    print("TEST 2: EC2 IMDSv1 metadata service (legacy)")
    print("=" * 60)
    try:
        import requests
        meta_resp = requests.get(
            'http://169.254.169.254/latest/meta-data/instance-id',
            timeout=2)
        meta_resp.raise_for_status()
        print(f"  Instance ID: {meta_resp.text}")
        print("  RESULT: DETECTED as AWS EC2 via IMDSv1")
        return True
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        print("  RESULT: NOT detected via IMDSv1")
        return False


def test_boto3_sts():
    """Test 3: boto3 STS GetCallerIdentity (works with IAM roles)."""
    print()
    print("=" * 60)
    print("TEST 3: boto3 STS GetCallerIdentity")
    print("=" * 60)
    try:
        import boto3
        import botocore.exceptions
        sts = boto3.client('sts', region_name='us-west-2')
        identity = sts.get_caller_identity()
        print(f"  Account: {identity['Account']}")
        print(f"  ARN:     {identity['Arn']}")
        print(f"  UserId:  {identity['UserId']}")
        print("  RESULT: DETECTED as AWS (valid AWS credentials)")
        return True
    except ImportError:
        print("  FAILED: boto3 not installed")
        print("  RESULT: Cannot test (install boto3)")
        return False
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        print("  RESULT: NOT detected via STS")
        return False


def test_boto3_region():
    """Test 4: boto3 session region / credential detection."""
    print()
    print("=" * 60)
    print("TEST 4: boto3 session metadata")
    print("=" * 60)
    try:
        import boto3
        import botocore.exceptions
        session = boto3.session.Session()
        region = session.region_name
        creds = session.get_credentials()
        print(f"  Region: {region}")
        if creds:
            resolved = creds.get_frozen_credentials()
            has_key = bool(resolved.access_key)
            has_token = bool(resolved.token)
            print(f"  Has access key: {has_key}")
            print(f"  Has session token: {has_token}")
            print(f"  Credential method: {creds.method}")
            if has_key:
                print("  RESULT: AWS credentials FOUND via boto3 session")
                return True
        print("  RESULT: No AWS credentials found")
        return False
    except ImportError:
        print("  FAILED: boto3 not installed")
        return False
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        return False


def test_env_vars():
    """Test 5: Check AWS-related environment variables."""
    print()
    print("=" * 60)
    print("TEST 5: AWS environment variables")
    print("=" * 60)
    aws_vars = [
        'AWS_DEFAULT_REGION',
        'AWS_REGION',
        'AWS_ACCESS_KEY_ID',
        'AWS_SECRET_ACCESS_KEY',
        'AWS_SESSION_TOKEN',
        'AWS_CONTAINER_CREDENTIALS_RELATIVE_URI',
        'AWS_CONTAINER_CREDENTIALS_FULL_URI',
        'AWS_EXECUTION_ENV',
        'AWS_LAMBDA_FUNCTION_NAME',
        'ECS_CONTAINER_METADATA_URI',
        'ECS_CONTAINER_METADATA_URI_V4',
    ]

    found_any = False
    for var in aws_vars:
        val = os.environ.get(var)
        if val:
            # Mask credentials
            if 'KEY' in var or 'SECRET' in var or 'TOKEN' in var:
                display = val[:8] + '...' if len(val) > 8 else '***'
            else:
                display = val
            print(f"  {var} = {display}")
            found_any = True
        else:
            print(f"  {var} = (not set)")

    if found_any:
        print("  RESULT: Some AWS env vars are set")
    else:
        print("  RESULT: No AWS env vars found")
    return found_any


def test_ecs_metadata():
    """Test 6: ECS/Fargate container metadata endpoint."""
    print()
    print("=" * 60)
    print("TEST 6: ECS/Fargate metadata endpoint")
    print("=" * 60)
    ecs_uri = os.environ.get('ECS_CONTAINER_METADATA_URI_V4') or \
              os.environ.get('ECS_CONTAINER_METADATA_URI')
    if not ecs_uri:
        print("  ECS metadata URI not in environment")
        print("  RESULT: NOT in ECS/Fargate container")
        return False
    try:
        import requests
        resp = requests.get(ecs_uri, timeout=2)
        resp.raise_for_status()
        data = resp.json()
        print(f"  Container: {data.get('DockerId', 'unknown')[:12]}")
        print("  RESULT: DETECTED as ECS/Fargate")
        return True
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        return False


def test_container_creds():
    """Test 7: AWS container credentials URI (ECS/EKS task role)."""
    print()
    print("=" * 60)
    print("TEST 7: Container credentials endpoint")
    print("=" * 60)
    rel_uri = os.environ.get('AWS_CONTAINER_CREDENTIALS_RELATIVE_URI')
    full_uri = os.environ.get('AWS_CONTAINER_CREDENTIALS_FULL_URI')
    if rel_uri:
        print(f"  Relative URI: {rel_uri}")
        print("  RESULT: Container role credentials available")
        return True
    elif full_uri:
        print(f"  Full URI: {full_uri[:40]}...")
        print("  RESULT: Container role credentials available")
        return True
    else:
        print("  No container credential env vars set")
        print("  RESULT: No container credentials endpoint")
        return False


def main():
    print("AWS Detection Diagnostic")
    print("Running from:", sys.executable)
    print("Platform:", sys.platform)
    print()

    results = {}
    results['IMDSv2'] = test_imdsv2()
    results['IMDSv1'] = test_imdsv1()
    results['boto3_STS'] = test_boto3_sts()
    results['boto3_session'] = test_boto3_region()
    results['env_vars'] = test_env_vars()
    results['ECS_metadata'] = test_ecs_metadata()
    results['container_creds'] = test_container_creds()

    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL"
        print(f"  {name:20s}: {status}")

    detected = any(results.values())
    print()
    if detected:
        # Which methods worked
        working = [k for k, v in results.items() if v]
        print(f"  --> AWS DETECTED via: {', '.join(working)}")
    else:
        print("  --> NOT detected as AWS by any method")

    return detected


if __name__ == '__main__':
    main()
