/*--------------------------------------------------------------------------*/
/*------------------------- File MCFSimplex.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Linear and Quadratic Min Cost Flow problems solver based on the (primal and
 * dual) simplex algorithm. Conforms to the standard MCF interface defined in
 * MCFClass.h.
 *
 * \Version 1.00
 *
 * \date 29 - 08 - 2011
 *
 * \author Alessandro Bertolini \n
 *         Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2008 - 2011 by Alessandro Bertolini, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFSimplex
 #define __MCFSimplex  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"
//Johannes Sommer

/*@} -----------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFSimplex_MACROS Compile-time switches in MCFSimplex.h
    There is only one macro in MCFSimplex, but it is very important!
    @{ */

#define QUADRATICCOST 0

/**< Setting QUADRATICCOST == 1 the solver can solve problems with linear and 
   quadratic costs too (but the latter only with the Primal Simplex).
   The reason for having a macro is that when quadratic costs are present the
   "arcType" struct has the additional field "quadraticCost" to hold it.
   Furthermore, the field "ident" is not created because the solver doesn't
   use the classical TLU tripartition. Instead, closed arcs and deleted arcs
   are characterized as follows:
   - closed arcs have the field "cost" to INFINITY (Inf<FNumber>());
   - deleted arcs have the field "upper" to INFINITY and the "tail" and
     "head" field are NULL.
   Furthermore, the solver needs the variables "ignoredEnteringArc" and
   "firstIgnoredEnteringArc", used to avoid nasty loops during the execution
   of the Quadratic Primal Simplex algorithm.
   If, instead, QUADRATICCOST == 0 then the solver can solve only problems
   with linear costs. Hence, the field "quadraticCost" is useless and it
   isn't created. Furthermore, Primal Simplex and Dual Simplex use the
   tripartition TLU to divide the arcs, so the solver creates the field
   "ident", which differentiates the set of the arcs in among deleted arcs,
   closed arcs, arcs in T, arcs in L, arcs in U.
   Thus, with QUADRATICCOST == 0 the solver cannot solve problems with
   quadratic costs, but it does solve problems with linear costs faster. */

/*@}  end( group( MCFCLASS_MACROS ) ) */ 

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFSimplex_CLASSES Classes in MCFSimplex.h
    @{ */

/** The MCFSimplex class derives from the abstract base class MCFClass, thus
    sharing its (standard) interface, and implements both the Primal and
    Dual network simplex algorithms for solving (Linear and Quadratic) 
    Min Cost Flow problems */

class MCFSimplex: public MCFClass 
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/** Public enum describing the parameters of MCFSimplex. */

 enum SimplexParam 
 { 
  kAlgPrimal = kLastParam , ///< parameter to set algorithm (Primal/Dual):
  kAlgPricing ,             ///< parameter to set algorithm of pricing
  kNumCandList ,            /**< parameter to set the number of candidate
			         list for Candidate List Pivot method */
  kHotListSize ,            /**< parameter to set the size of Hot List
			         for Candidate List Pivot method */
  kRecomputeFOLimits ,      /**< parameter to set the number of iterations
                                 in which quadratic Primal Simplex computes 
                                 "manually" the f.o. value */
  kEpsOpt                   /**< parameter to set the precision of the 
                                 objective function value for the 
				 quadratic Primal Simplex */
  };
    
/** Public enum describing the pricing rules in MCFSimplex::SetAlg(). */

 enum enumPrcngRl 
 { 
  kDantzig = 0,        ///< Dantzig's rule (most violated constraint)
  kFirstEligibleArc ,  ///< First eligible arc in round-robin
  kCandidateListPivot  ///< Candidate List Pivot Rule
  };

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

 MCFSimplex( cIndex nmx = 0 , cIndex mmx = 0 );

/**< Constructor of the class, as in MCFClass::MCFClass(). */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
	       cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
	       cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
	       cIndex_Set pEn = NULL );

/**< Inputs a new network, as in MCFClass::LoadNet(). */

/*--------------------------------------------------------------------------*/

 void SetAlg( bool UsPrml , char WhchPrc );

/**< Selects which algorithm (Primal vs Dual Network Simplex), and which
   pricing rule within the algorithm, is used.

   If UsPrml == TRUE then the Primal Network Simplex algorithm is used,
   otherwise the Dual Network Simplex is used.

   The allowed values for WhchPrc are:

   - kDantzig            Dantzig's pricing rule, i.e., most violated dual 
                         constraint; this can only be used with the Primal 
                         Network Simplex

   - kFirstEligibleArcA  "dumb" rule, first eligible arc in round-robin;

   - kCandidateListPivot Candidate List Pivot Rule

   If this method is *not* called, the Primal Network Simplex with the
   Candidate List Pivot Rule (the best setting on most instances) is
   used. */

/*--------------------------------------------------------------------------*/

 void SetPar( int par , int val );

/**< Set general integer parameters.

   @param par   is the parameter to set; since this method accepts an int
                value, the enum SimplexParam can be used in addition to the
                enum MCFParam to specify the integer parameters (every enum
		is an int).

   @param value is the value to assign to the parameter. */

/*-------------------------------------------------------------------------*/

 void SetPar( int par , double val );

/**< Set general float parameters.

   @param par   is the parameter to set; since this method accepts an int
                value, the enum SimplexParam can be used in addition to the
                enum MCFParam to specify the float parameters (every enum
		is an int).

   @param value is the value to assign to the parameter. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

 void SolveMCF( void );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 void MCFGetX( FRow F , Index_Set nms = NULL ,
	       cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

 void MCFGetRC( CRow CR , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() ) ;

 CNumber MCFGetRC( cIndex i );

/*--------------------------------------------------------------------------*/

 void MCFGetPi( CRow P , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

 FONumber MCFGetFO( void );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

 virtual void MCFArcs( Index_Set Startv , Index_Set Endv ,
		       cIndex_Set nms = NULL , cIndex strt = 0 ,
		       Index stp = Inf<Index>() );

 inline Index MCFSNde( cIndex i );

 inline Index MCFENde( cIndex i );

/*--------------------------------------------------------------------------*/
  
 void MCFCosts( CRow Costv , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 inline CNumber MCFCost( cIndex i );

/*--------------------------------------------------------------------------*/

 void MCFQCoef( CRow Qv , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 inline CNumber MCFQCoef( cIndex i );

/*--------------------------------------------------------------------------*/

 void MCFUCaps( FRow UCapv , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 inline FNumber MCFUCap( cIndex i );

/*--------------------------------------------------------------------------*/

 void MCFDfcts( FRow Dfctv , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 inline FNumber MCFDfct( cIndex i );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*------- Changing the costs, QCoef, deficits and upper capacities ---------*/
/*--------------------------------------------------------------------------*/

 void ChgCosts( cCRow NCost , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 void ChgCost( Index arc , cCNumber NCost );

/*--------------------------------------------------------------------------*/

 void ChgQCoef( cCRow NQCoef = NULL , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 void ChgQCoef( Index arc , cCNumber NQCoef );

/*--------------------------------------------------------------------------*/

 void ChgDfcts( cFRow NDfct , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 void ChgDfct( Index nod , cFNumber NDfct );

/*--------------------------------------------------------------------------*/

 void ChgUCaps( cFRow NCap , cIndex_Set nms = NULL ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

 void ChgUCap( Index arc , cFNumber NCap );

/*--------------------------------------------------------------------------*/
/*--------------- Modifying the structure of the graph ---------------------*/
/*--------------------------------------------------------------------------*/

 void CloseArc( cIndex name );

 void DelNode( cIndex name );

 bool IsClosedArc( cIndex name );

 void OpenArc( cIndex name ) ;

 Index AddNode( cFNumber aDfct );

 void ChangeArc( cIndex name , cIndex nSS = Inf<Index>() ,
		 cIndex nEN = Inf<Index>() );

 void DelArc( cIndex name );

 Index AddArc( cIndex Start , cIndex End , cFNumber aU , cCNumber aC );

 bool IsDeletedArc( cIndex name );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

 ~MCFSimplex();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE DATA TYPES ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Let T \subseteq A be a spanning tree, and consider some node v \in V --*/
/*-- \setminus { 0 }. There is an unique (undirected) path, denoted by    --*/
/*-- P(v), defined by T from v to the root node 0. The arc in P(v), which --*/
/*-- is incident to v, is called the *basic arc* of v. The other terminal --*/
/*-- node u of this basic arc is called the *father node* of v. The       --*/
/*-- spanning tree T is represented saving the basic arc of every node,   --*/
/*-- and maintaining the order of the nodes and the depth as to T root    --*/
/*-- after a Post-Visit of T. This order is saved in a bidirectional list --*/
/*-- written in the node.                                                 --*/
/*--                                                                      --*/
/*-- The Primal Simplex uses a different data structure than the Dual     --*/
/*-- Simplex, because the Dual Simplex needs additional data (mainly the  --*/
/*-- Backward Star and Forward Star. Furthermore, the Primal Simplex uses --*/
/*-- different data structures in the quadratic case.                     --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 struct arcPType;   // pre-declaration of the arc structure (pointers to arcs
                    // are contained in the node structure) for Primal Simplex

 struct arcDType;   // pre-declaration of the arc structure (pointers to arcs
                    // are contained in the node structure) for Dual Simplex

 typedef double iteratorType;  // type for the iteration counter and array
                               // "whenInT2"

 struct nodePType {       // node structure for Primal Simplex - - - - - - - -
  nodePType *prevInT;     // previous node in the order of the Post-Visit on T

  nodePType *nextInT;     // next node in the order of the Post-Visit on T

  arcPType *enteringTArc; // entering basic arc of this node     

  FNumber balance;        // supply/demand of this node; a node is called a
                          // supply node, a demand node, or a transshipment
                          // node depending upon whether balance is larger
                          // than, smaller than, or equal to zero
  #if( QUADRATICCOST )
   CNumber sumQuadratic;  // the sum of the quadratic coefficients of the tree's arcs
                          // from root of T to the node

   FONumber potential;    // the node potential corresponding with the flow
                          // conservation constrait of this node
  #else
   CNumber potential;      // the node potential corresponding with the flow
                           // conservation constrait of this node
  #endif

  int subTreeLevel;        // the depth of the node in T as to T root

  };                       // end( struct( nodePType ) )

 struct nodeDType {       // node structure for Dual Simplex - - - - - - - - -
  nodeDType *prevInT;     // previous node in the order of the Post-Visit on T

  nodeDType *nextInT;     // next node in the order of the Post-Visit on T

  arcDType *enteringTArc; // entering basic arc of this node     

  FNumber balance;        // supply/demand of this node; a node is called a
                          // supply node, a demand node, or a transshipment
                          // node depending upon whether balance is larger
                          // than, smaller than, or equal to zero

 #if( QUADRATICCOST )
  CNumber sumQuadratic;   // the sum of the quadratic coefficients of the tree's arcs
                          // from root of T to the node
        
  FONumber potential;     // the node potential corresponding with the flow
                          // conservation constrait of this node
 #else
  CNumber potential;      // the node potential corresponding with the flow
                          // conservation constrait of this node
 #endif

  int subTreeLevel;       // the depth of the node in T as to T root
  iteratorType whenInT2;  // the last iteration where a node is in subtree T2

  Index numArcs;          // the number of the arcs which enter/exit from node
  arcDType *firstBs;      // the first arc in the node's Backward Star
  arcDType *firstFs;      // the first arc in the node's Forward Star

  };                      // end( struct( nodeDType ) )

 struct arcPType {        // arc structure for Primal Simplex - - - - - - - -
  nodePType *tail;        // tail node
  nodePType *head;        // head node
  FNumber flow;           // arc flow
  CNumber cost;           // arc linear cost

  #if( QUADRATICCOST )
   CNumber quadraticCost; // arc quadratic cost
  #else
   char ident;            // tells if arc is deleted, closed, in T, L, or U
  #endif

  FNumber upper;          // arc upper bound
  };                      // end( struct( arcPType ) )

 struct arcDType {        // arc structure for Primal Simplex - - - - - - - -
  nodeDType *tail;        // tail node
  nodeDType *head;        // head node
  FNumber flow;           // arc flow
  CNumber cost;           // arc linear cost

  #if( QUADRATICCOST )
   CNumber quadraticCost; // arc quadratic cost
  #else
   char ident;            // indicates if arc is deleted, closed, in T, in L, or in U
  #endif

  FNumber upper;          // arc upper bound
  arcDType *nextBs;       // the next arc in the Backward Star of the arc's head
  arcDType *nextFs;       // the next arc in the Forward Star of the arc's tail

  };                      // end( struct( arcDType ) )

 struct primalCandidType {  // Primal Candidate List- - - - - - - - - - - - -
  arcPType *arc;            // pointer to the violating primal bound arc

  #if(QUADRATICCOST )
   FONumber absRC;          // absolute value of the arc's reduced cost
  #else
   CNumber absRC;           // absolute value of the arc's reduced cost
  #endif
  };                        // end( struct( primalCandidateType ) )

 struct dualCandidType {    // Dual Candidate List- - - - - - - - - - - - - -
  nodeDType *node;          //  deepest node violating the dual bound arc
  FNumber absInfeas;        // absolute value of the arc's flow infeasibility
  };                        // end( struct( dualCandidateType ) )

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 bool usePrimalSimplex;         // TRUE if the Primal Network Simplex is used

 char pricingRule;              // which pricing rule is used

 nodePType *nodesP;             // vector of nodes: points to the n + 1 node
                                // structs (including the dummy root node)
                                // where the first node is indexed by zero
                                // and the last node is the dummy root node

 nodePType *dummyRootP;         // the dummy root node

 nodePType *stopNodesP;         // first infeasible node address = nodes + n 
  
 arcPType *arcsP;               // vector of arcs; this variable points to
                                // the m arc structs.

 arcPType *dummyArcsP;          // vector of artificial dummy arcs: points
                                // to the artificial dummy arc variables and
                                // contains n arc structs. The artificial
                                // arcs are used to build artificial feasible
                                // starting bases. For each node i, there is
                                // exactly one dummy arc defined to connect
                                // the node i with the dummy root node.
    
 arcPType *stopArcsP;           // first infeasible arc address = arcs + m 

 arcPType *stopDummyP;          // first infeasible dummy arc address
                                // = arcs + m + n

 arcPType *arcToStartP;         // Dantzig Rule and First Eligible Arc Rule
                                // start their search from this arc

 nodeDType *nodesD;             // vector of nodes: points to the n + 1 node
                                // structs (including the dummy root node)
                                // where the first node is indexed by zero
                                // and the last node is the dummy root node

 nodeDType *dummyRootD;         // the dummy root node

 nodeDType *stopNodesD;         // first infeasible node address = nodes + n 
  
 arcDType *arcsD;               // vector of arcs; this variable points to
                                // the m arc structs.

 arcDType *dummyArcsD;          // vector of artificial dummy arcs: points to
                                // to the artificial dummy arc variables and
                                // contains n arc structs. The artificial
                                // arcs are used to build artificial feasible
                                // starting bases. For each node i, there is
                                // exactly one dummy arc defined to connect
                                // the node i with the dummy root node.

 arcDType *stopArcsD;           // first infeasible arc address = arcs + m 

 arcDType *stopDummyD;          // first infeasible dummy arc address
                                // = arcs + m + n

 arcDType *arcToStartD;         // Dantzig Rule and First Eligible Arc Rule
                                // start their search from this arc

 iteratorType iterator;         // the current number of iterations
    
 primalCandidType *candP;       // every element points to an element of the
                                // arcs vector which contains an arc violating 
                                // dual bound

 dualCandidType *candD;         // every element points to an element of the
                                // arcs vector which contains an arc violating 
                                // primal bond

 Index numGroup;                // number of the candidate lists
    
 Index tempCandidateListSize;   // hot list dimension (it is variable)
    
 Index groupPos;                // contains the actual candidate list
    
 Index numCandidateList;        // number of candidate lists
    
 Index hotListSize;             // number of candidate lists and hot list dimension

 Index forcedNumCandidateList;  // value used to force the number of candidate list
    
 Index forcedHotListSize;       // value used to force the number of candidate list
                                // and hot list dimension

 bool newSession;               // true if algorithm is just started
  
 CNumber MAX_ART_COST;          // large cost for artificial arcs

 FNumber *modifiedBalance;      // vector of balance used by the PostVisit

 FONumber EpsOpt;               // the precision of the objective function value
                                // for the quadratic case of the Primal Simplex

 int recomputeFOLimits;         // after how many iterations the quadratic Primal
                                // Simplex computes "manually" the o.f. value

 #if( QUADRATICCOST )
  FONumber foValue;             // the temporary objective function value
 #endif

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void MemAlloc( void );

/**< Method to allocate memory for the main data structures. It creates the
   vector of nodes, the vector of arcs (real and dummy). If the algorithm is
   using the Dual Simplex, it creates also the vectors whenInT2, firstIn and
   nextIn, usefull to identify the next entering. */

/*--------------------------------------------------------------------------*/

  void MemDeAlloc( bool whatDeAlloc );

/**< Method to deallocate memory for the main data structures created in
   MemAlloc(). */

/*--------------------------------------------------------------------------*/

  void MemAllocCandidateList( void );

/**< Method to allocate memory for the data structure used by the Candidate
   List Pivot Rule. It creates the vector of candP (or candD in the Dual
   Simplex case), determining the size of this vector on the basis of number
   of arcs or nodes. */

/*--------------------------------------------------------------------------*/

  void MemDeAllocCandidateList( void );

/**< Method to deallocate memory for the data structures created in
   MemAllocCandidateList(). */

/*--------------------------------------------------------------------------*/

  void CreateInitialPrimalBase( void );

/**< Method to create an initial feasible primal base. Add one node (dummyRoot)
   to the network and connect this dummy root to each real node with n dummy
   arcs. For each source node i, a dummy arc (i, r) exists. For each sink or
   transit node j, a dummy arc (r, j) exists. The source nodes send their flows
   to the dummy root through their dummy arcs, so the dummy root send all to
   the sink nodes. The dummy root balance is 0, the costs of dummy arcs are
   fixed to "infinity". */

/*--------------------------------------------------------------------------*/

  void CreateInitialDualBase( void );

/**< Method to create an initial feasible dual base. Add one node (dummyRoot)
   to the network and connect this dummy root to each real node with n dummy
   arcs. The source nodes send their flows to the dummy root through their
   dummy arcs, so the dummy root send all to the sink nodes. The dummy root
   balance is 0, the costs of dummy arcs are fixed to "infinity". */

/*--------------------------------------------------------------------------*/

  void CreateAdditionalDualStructures( void );

/**< The Dual Simplex needs nodes' Backward and Forward Star to work. So when
   the Primal Simplex runs, these structures don't exist. When the Dual
   Simplex starts, these structure are created in this method. */

/*--------------------------------------------------------------------------*/

  void PrimalSimplex( void );

/**< Main method to implement the Primal Simplex algorithm. */

/*--------------------------------------------------------------------------*/

  void DualSimplex( void );

/**< Main method to implement the Dual Simplex algorithm. */

/*--------------------------------------------------------------------------*/

  template<class N, class A>
  void UpdateT( A *h , A *k , N *h1 , N *h2 , N *k1 , N *k2 );

  /**< Method to update the spanning tree T every iteration.
     The spanning tree T is implemented with a bidirectional list stored
     in the node structure, which represents the nodes' order after a
     Post-Visit. The parameter of the method are the outgoing arc "h", the
     incoming arc "k", and the four nodes on the extremity of these arcs
     (for example h2 is the deepest node of the outgoing arc). Removing the
     arc "h" splits T in two subtrees: T1 (which contains the root of T) and
     T2, which will be re-connected by the incoming arc "k".
     T2 will be reordered; in fact, the node "k2" becomes the root of T2 
     instead of "h2" and the hierarchy of T2 will be overturned. Then T2 will
     be moved; the root of T2 is changed, therefore the predecessor of the
     root will become the node "k1".

     This method uses the methods cutAndUpdateSubtree() and pasteSubtree().
     First it cuts the node "k2" and its subtree from T using
     cutAndUpdateSubtree(). Moreover the method cutAndUpdateSubtree()
     updates the field "subTreeLevel" of every subtree's nodes, since k2's
     subtree will be moved from the bottom to the top of T2. Then the method
     pasteSubtree() puts this subtree in the bidirectional list after the node
     "k1". The same operations will be applied to the old precedessor of "k2" 
     (which will become one of the childs of "k2"). This second subtree will
     be cut, the subTreeLevel fields will be updated, and it will be inserted
     in the bidirectional list after the k2's subtree. This is iterated until
     the node "h2" is reached. */

/*--------------------------------------------------------------------------*/

  template<class N>
  N* CutAndUpdateSubtree( N *root, int delta );

/**< This method cuts a generic subtree from the spanning tree T. Then it
   updates the field "subTreeLevel" of every subtree's nodes adding the value
   "delta". This method returns the last node of the subtree. */

/*--------------------------------------------------------------------------*/

  template<class N>
  void PasteSubtree( N *root , N *lastNode , N *previousNode );

/**< This method inserts a generic subtree with root passed by parameter into
   the spanning tree T, between the nodes "previousNode" and "lastNode". */

/*--------------------------------------------------------------------------*/

  arcPType* RuleDantzig( void );

/**< This method returns an arc which violates the dual conditions. It searchs
   the arc with most violation of dual conditions in the entire set of real
   arcs. It can be used only by the Primal Simplex in the case of networks
   with linear costs. */

/*--------------------------------------------------------------------------*/

  arcPType* PRuleFirstEligibleArc( void );

/**< This method returns the first found arc which violates the dual conditions
   in the case of Primal Simplex, the primal condition in the case of Dual
   Simplex. It can be used only in the case of networks with linear costs. */

/*--------------------------------------------------------------------------*/

  arcDType* DRuleFirstEligibleArc( void );

/**< This method returns the first found arc which violates the dual conditions
   in the case of Primal Simplex, the primal condition in the case of Dual
   Simplex. It can be used only in the case of networks with linear costs. */

/*--------------------------------------------------------------------------*/

  arcPType* RulePrimalCandidateListPivot( void );

/**< This method returns an arc which violates the dual conditions. It searches
   the arc with most violation of dual conditions in a small set of candidate
   arcs. In every iteration the method rebuilds this set of arcs executing three
   phases:

   - in the first phase it analyzes the remaining arcs and delete the arcs which
     don't violate the dual condition any more;

   - in the second phase it tries to fill the set, so it searchs other arcs 
     which violate the dual condition: the set of arcs is divided into "buckets"
     which are searched sequentially until the candidate list is full; the
     last visited bucket is retained, and the search is restarted from that
     one at later iterations

   - in the third phase the small set of candidate arcs is ordered according
     to the violation of dual condition by the method SortPrimalCandidateList() 
     using an implementation of the algorithm "quicksort".

    At last the method returns the first arc in the ordered small set. If the
    arc doesn't exist (the set is empty), it returns NULL. */

/*--------------------------------------------------------------------------*/

  inline void InitializePrimalCandidateList( void );

/**< Method to initialize some important details for Primal Candidate List
   Rule. */

/*--------------------------------------------------------------------------*/

  inline void SortPrimalCandidateList( Index min , Index max );

/**< Method to order the little set of candidate arcs according to
   infeasibility of dual conditions. It implements the "quicksort"
   algorithm.  */

/*--------------------------------------------------------------------------*/

  arcDType* RuleDualCandidateListPivot( void );

/**< Similar to RulePrimalCandidateListPivot() for the Dual Simplex. */

/*--------------------------------------------------------------------------*/

  inline void InitializeDualCandidateList( void );

/**< Method to initialize some important details for Dual Candidate List Rule.
    */

/*--------------------------------------------------------------------------*/

  inline void SortDualCandidateList( Index min , Index max );

/**< Similar to SortPrimalCandidateList() for the Dual Simplex. */

/*--------------------------------------------------------------------------*/

  template<class N, class RCT>
  inline void AddPotential( N *r , RCT delta );

/**< Method to quickly update the dual solutions. During the change of the
   base, the potential of node "k2" (deepest node in T of incoming arc "k",
   and new root of T2) changes according to the new structure of T. In fact,
   the precedessor of "k2" changes: now the predecessor of "k2"is "k1" (the
   other node of the incoming arc "k"). This change of predecessor causes a
   change of potential of "k2". The change of potential of "k2" causes the
   changes of potential of all nodes of T2. This method computes the change
   of potential of "k2" and applies it on all the nodes of T2. */

/*--------------------------------------------------------------------------*/

  template<class N>
  inline void ComputePotential( N *r );

/**< Method to update the dual solutions. It computes all the potential of the
   nodes of the subtree which has r as root. */

/*--------------------------------------------------------------------------*/

  inline void ResetWhenInT2( void );

/**< Method to order the small set of candidate arcs according to dual
   infeasibility. It implements the algorithm "quicksort". */

/*--------------------------------------------------------------------------*/

  void CreateInitialPModifiedBalanceVector( void );

/**< Method to initialize the vector of modified balance for the postvisit on
   T in the Primal Simplex data structure. */

/*--------------------------------------------------------------------------*/
    
  void PostPVisit( nodePType *r );

/**< Method to calculate the flow on the basic arcs with the Primal Simplex's
   data structure. It uses the set of the upper bound arcs, the construction
   of a modified balance vector and the postvisit on T. */

/*--------------------------------------------------------------------------*/

  void BalanceFlow( nodePType *r );

/**< This method works after the method PostPVisit (in a Primal context). It
   restores primal admissimibility on the r's subtree. It starts from the leaf
   of the subtree and goes up to the root, using the method AdjustFlow. */

/*--------------------------------------------------------------------------*/

  void AdjustFlow( nodePType *r );

/**< This method checks the primal admissibility of the r's basic entering arc.
    If it is out of bounds, the method removes it from the tree (and keeps the
    relative dummy arc) and push flow in the cycle (some tree's arc and the old
    entering arc) to restores the right balances of the node. */

/*--------------------------------------------------------------------------*/

  void CreateInitialDModifiedBalanceVector( void );

/**< Method to initialize the vector of modified balance for the postvisit on
   T in the Dual Simplex data structure. */

/*--------------------------------------------------------------------------*/
    
  void PostDVisit( nodeDType *r );

/**< Method to calculate the flow on the basic arcs with the Dual Simplex's data
   structure, using the set of the upper bound arcs, the construction of a
   modified balance vector and the postvisit on T. */

/*--------------------------------------------------------------------------*/

  template<class N, class A>
  inline N* Father( N *n, A *a );

/**< Method to find the predecessor of the node in the tree. */

/*--------------------------------------------------------------------------*/

 #if(QUADRATICCOST)
  template<class A>
  inline FONumber ReductCost( A *a );
 #else
  template<class A>
  inline CNumber ReductCost( A *a );
 #endif

/**< Method to calculate the reduct cost of the arc. */

/*--------------------------------------------------------------------------*/

  inline FONumber GetFO( void );

/**< Method to calculate the temporary (or the final) objective function
   value. */

/*--------------------------------------------------------------------------*/
    
  void PrintPNode( nodePType *nodo );

/**< Method to print the "name" of the node in the Primal Simplex. */

/*--------------------------------------------------------------------------*/

  void PrintPArc( arcPType *arc );

/**< Method to print the "name" of the arc in the Primal Simplex. */

/*--------------------------------------------------------------------------*/

  void PrintDNode( nodeDType *nodo );

/**< Method to print the "name" of the node in the Dual Simplex. */

/*--------------------------------------------------------------------------*/

  void PrintDArc( arcDType *arc );

/**< Method to print the "name" of the arc in the Dual Simplex. */

/*--------------------------------------------------------------------------*/

  nodePType* RecoverPNode( Index ind );

/**< Method to find a node (in the Primal Simplex) using its index. */

/*--------------------------------------------------------------------------*/

  arcPType* RecoverPArc( nodePType *tail , nodePType *head );

/**< Method to find an arc (in the Primal Simplex) using 2 pointers to tail
   node and head node. */

/*--------------------------------------------------------------------------*/

  nodeDType* RecoverDNode( Index ind );

/**< Method to find a node (in the Dual Simplex) using its index. */

/*--------------------------------------------------------------------------*/

  arcDType* RecoverDArc( nodeDType *tail , nodeDType *head ); 

/**< Method to find an arc (in the Dual Simplex) using 2 pointers to tail
   node and head node. */

/*--------------------------------------------------------------------------*/

  void infoPNode( nodePType *node , int tab );

/**< Method to print some information of the node (in the Primal Simplex). */

/*--------------------------------------------------------------------------*/

  void infoPArc( arcPType *arc , int ind , int tab ); 

/**< Method to print some information of the arc (in the Primal Simplex). */

/*--------------------------------------------------------------------------*/

  void infoDNode( nodeDType *node , int tab );

/**< Method to print some information of the node (in the Dual Simplex). */

/*--------------------------------------------------------------------------*/

  void infoDArc( arcDType *arc , int ind , int tab );

/**< Method to print some information of the arc (in the Dual Simplex). */

/*--------------------------------------------------------------------------*/

  void ShowSituation( int tab );

/**< Method to show the actual complete situation. */

/*--------------------------------------------------------------------------*/
    
  };  // end( class MCFSimplex )

/* @} end( group( MCFSimplex_CLASSES ) ) */

#endif  /* MCFSimplex.h included */

/*-------------------------------------------------------------------------*/
/*-------------------inline methods implementation-------------------------*/
/*-------------------------------------------------------------------------*/

inline MCFClass::Index MCFSimplex::MCFSNde( MCFClass::cIndex i )
{
 if( usePrimalSimplex )
  return( Index( ( (arcsP + i)->tail - nodesP + 1 ) - USENAME0 ) );
 else
  return( Index( ( (arcsD + i)->tail - nodesD + 1 ) - USENAME0 ) );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::Index MCFSimplex::MCFENde( MCFClass::cIndex i )
{
 if( usePrimalSimplex )
  return( Index( ( (arcsP + i)->head - nodesP + 1 ) - USENAME0 ) );
 else
  return( Index( ( (arcsD + i)->head - nodesD + 1 ) - USENAME0 ) );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::CNumber MCFSimplex::MCFCost( MCFClass::cIndex i )
{
 if( usePrimalSimplex )
  return( (arcsP + i)->cost );
 else
  return( (arcsD + i)->cost );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::CNumber MCFSimplex::MCFQCoef( MCFClass::cIndex i )
{
 #if( QUADRATICCOST )
  if( usePrimalSimplex )
   return( (arcsP + i)->quadraticCost );
  else
   return( (arcsD + i)->quadraticCost );
 #else
  return( 0 );
 #endif
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::FNumber MCFSimplex::MCFUCap( MCFClass::cIndex i )
{
 if( usePrimalSimplex )
  return( (arcsP + i)->upper );
 else
  return( (arcsD + i)->upper );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::FNumber MCFSimplex::MCFDfct( MCFClass::cIndex i )
{
 if( usePrimalSimplex )
  return( (nodesP + i)->balance );
 else
  return( (nodesD + i)->balance );
 }

/*-------------------------------------------------------------------------*/

#if( QUADRATICCOST )

template <class A>
inline MCFSimplex::FONumber MCFSimplex::ReductCost( A *a )
{
 FONumber redc = (a->tail)->potential - (a->head)->potential;
 redc = redc + a->cost;
 redc = redc + a->quadraticCost * a->flow;
 return( redc );
 }

#else

template <class A>
inline MCFSimplex::CNumber MCFSimplex::ReductCost( A *a )
{
 CNumber redc = (a->tail)->potential - (a->head)->potential;
 redc = redc + a->cost;
 return( redc );
 }

#endif

/*-------------------------------------------------------------------------*/
 
#if( OPT_USE_NAMESPACES )
};  // end( namespace MCFClass_di_unipi_it )
#endif

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/*---------------------- End File MCFSimplex.h ----------------------------*/
/*-------------------------------------------------------------------------*/
