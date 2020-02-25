/*--------------------------------------------------------------------------*/
/*----------------------------- File SPTree.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of several "classic" Shortest Path Tree algorithms to
 * solve uncapacitated single-source Min Cost Flow problems. The actual
 * algorithm can be chosen at compile time by setting a proper switch.
 * Conforms to the standard MCF interface defined in MCFClass.h.
 *
 * \version 1.96
 *
 * \date 16 - 10 - 2018
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 1996 - 2018 by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _SPTree
 #define _SPTree  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*---------------------------- MACROS --------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SPTREE_MACROS Compile-time switches in SPTree.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

/*------------------------------ SPT_ALGRTM --------------------------------*/

#define SPT_ALGRTM 0

/**< This macro decides which SPT algorithm has to be used.
   Possible values are:

   - 0  =>  LQueue
   - 1  =>  LDeque
   - 2  =>  (currently unused)
   - 3  =>  Dijkstra
   - 4  =>  Heap

   For algorithms based on priority lists, the macro LABEL_SETTING [see
   below] can be set to 1 to say that the algorithm is of the
   "label-setting" (nodes only exit from Q once) rather than of the
   "label-correcting" (nodes may exit from Q more than once) type. */

#if( SPT_ALGRTM <= 2 )
 #define LABEL_SETTING 0
 ///< this is a label-correcting SPT algorithm
#else
 #define LABEL_SETTING 1

 /**< This macro decides if the "label-setting" style is used.

    With a priority lists, the SPT algorithm applied to SPT problems
    with *all nonnegative arc costs* has the "label-setting" property:
    nodes only exit from Q once, hence when a node exits from Q its
    label is permanently set.

    If LABEL_SETTING > 0 the code will assume that this property holds
    and implement some things accordingly; in particular, the algorithm
    is terminated when the last destination is extracted from Q even
    though Q is still nonempty.

    \warning Solving a SPT algorithm with negative arc costs with
             LABEL_SETTING > 0 may produce a suboptimal solution. */

 #if( SPT_ALGRTM == 4 )
  #define HeapCard 2

  /**< Number of sons of each node in the heap.
     SPT_ALGRTM == 4 means using a C-ary heap to hold the node set Q: each
     HeapCard is the ariety of the heap, i.e. the max number of sons of a
     node in the heap. Special treatment is deserved to the case
     HeapCard == 2. */
 #endif
#endif

/*------------------------------ SPT_STRTN ---------------------------------*/

#define SPT_STRTN 1

/* Decides if the "start node" information for each arc is explicitly kept
   in a data structure. If SPT_STRTN == 0 the "start node" information is
   computed in O( ln( n ) ) each time it is needed. This is only used in
   methods for reading or changing the data of the problem, and not in the
   main (SPT) algorithm, so it may not be too costly. If SPT_STRTN == 1
   instead the data structure is constructed; if SAME_GRPH_SPT > 0, the
   data structure is "static".

   \note This switch does not appear in the manual because the current
         implementation of Startn() for SPT_STRTN == 0 is flawed. */

/*------------------------------ ORDRD_NMS ---------------------------------*/

#define ORDRD_NMS 1

/**< Decides if arc names in MCFGetX() are ordered.
   If ORDRD_NMS > 0, and MCFGetX() [see below] is asked for a "sparse" flow
   solution (i.e., nms != 0), then the set of indices returned at the end
   of the method is ordered in increasing sense. If ORDRD_NMS == 0 instead,
   the set of indices may not be ordered.

   ORDRD_NMS > 0 may be useful for some applications, but it is more costly
   (basically, it requires either to compute the "dense" flow solution or to
   sort a vector). Also, "sparse" flow solutions in this class are guaranteed
   to contain no more than n - 1 nonzeroes, hence if ORDRD_NMS == 0 then the
   parameter `F' in MCFGetX( F , nms ) can actually point to a (n - 1)-vector,
   while if ORDRD_NMS > 0 it must point to a m-vector anyway. */

/*----------------------------- SAME_GRPH_SPT ------------------------------*/

#define SAME_GRPH_SPT 0

/**< Decides if all MCFClass instances share the same graph.
   If SAME_GRPH_SPT > 0, then all the instances of the class will work on the
   same "topological" network, while the costs, capacities and supplies can
   change from one instance to another. This allows implementations to share
   some data structures describing the graph, e.g. by declaring them "static",
   saving memory when multiple instances of the solver are active at the same
   time. */

/*----------------------------- DYNMC_MCF_SPT ------------------------------*/

#define DYNMC_MCF_SPT 0

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_SPT > 0, some of the methods of the public interface of
   class that allow to change the topology of the underlying network are
   actually implemented. Possible values of this macro are:

   - 0 => the topology of the graph cannot be changed;

   - 1 => all the methods that change the topology of the graph are
          implemented. */

/*@} -----------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SPTREE_CLASSES Classes in SPTree.h
    @{ */

/** The SPTree class derives from the abstract base class MCFClass, thus
    sharing its (standard) interface, and implements Shortest Path Tree
    algorithms for solving "uncapacitated" (Linear) Min Cost Flow
    problems with one source node.

    \warning The SPT algorithm will enter in an infinite loop if a directed
             cycle of negative cost exists in the graph: there is no check
	     about this in the code. */

class SPTree : public MCFClass
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
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   SPTree( cIndex nmx = 0 , cIndex mmx = 0 , bool Drctd = true );

/**< Constructor of the class.

   For the meaning of nmx and mmx see MCFClass::MCFClass().

   The parameter `Drctd' tells if the given graph has really to be
   understood as directed (default), i.e., if the i-th arc is
   Sn[ i ] --> En[ i ], or undirected, i.e., the i-th arc is
   Sn[ i ] <--> En[ i ]. Undirected graphs are internally implemented by
   doubling each arc, but this is completely hidden by the interface. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = 0 , cCRow pC = 0 ,
		 cFRow pDfct = 0 , cIndex_Set pSn = 0 ,
		 cIndex_Set pEn = 0 ) override;

/**< Inputs a new network, as in MCFClass::LoadNet().

   Arcs with pC[ i ] == Inf<CNumber>() do not "exist". If DYNMC_MCF_SPT > 0,
   these arcs are "closed".

   If DYNMC_MCF_SPT == 0 but SAME_GRPH_SPT > 0, these arcs are dealt with
   explicitly, and can be put back into the formulation by simply changing
   their cost. Note that, however, this is less efficient than eliminating
   them explicitly from the problem.

   If DYNMC_MCF_SPT == 0 and SAME_GRPH_SPT == 0, these arcs are just removed
   from the formulation. However, they have some sort of a "special status"
   (after all, if the user wants to remove them completely he/she can just
   change the data), in that they are still counted into the number of arcs
   of the graph and they will always have 0 flow and Inf<CNumber>() reduced
   cost as "closed" or "deleted" arcs. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   using MCFClass::MCFGetX;  // the ( void ) method, which is otherwise hidden

   void MCFGetX( FRow F , Index_Set nms = 0  ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*--------------------------------------------------------------------------*/

   using MCFClass::MCFGetRC;  // the ( void ) method, which is otherwise hidden

   void MCFGetRC( CRow CR , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline CNumber MCFGetRC( cIndex i ) override;

/*--------------------------------------------------------------------------*/

   void MCFGetPi( CRow P , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   cCRow MCFGetPi( void ) override;

/**< Same meaning as MCFClass::MCFGetPi().

   \note Some of the potentials may be + Inf<CNumber>(): this means that

   - the node is *not* a destination and it cannot be reached from the Origin
     (however, this does *not* mean that the problem is unfeasible);

   - if LABEL_SETTING == 1, the node is *not* a destination and it has not
     been reached during the algorithm. */

/*--------------------------------------------------------------------------*/

   SPTree::FONumber MCFGetFO( void ) override;

/**< Same meaning as MCFClass::MCFGetFO().

   \note if not all the specified destinations can be reached from the
         Origin, returns Inf<FONumber>(). */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = 0  ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline Index MCFSNde( cIndex i ) override;

   inline Index MCFENde( cIndex i ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void MCFCosts( CRow Costv , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline CNumber MCFCost( cIndex i ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void MCFUCaps( FRow UCapv , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline FNumber MCFUCap( cIndex i ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void MCFDfcts( FRow Dfctv , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline FNumber MCFDfct( cIndex i ) override;

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgCost( Index arc , cCNumber NCost ) override;

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgDfct( Index nod , cFNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap , cIndex_Set nms = 0  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgUCap( Index arc , cFNumber NCap ) override;

/*--------------------------------------------------------------------------*/
/*--------------- Modifying the structure of the graph ---------------------*/
/*--------------------------------------------------------------------------*/

  void CloseArc( cIndex name ) override;

  inline bool IsClosedArc( cIndex name ) override;

  void DelNode( cIndex name ) override;

  void OpenArc( cIndex name ) override;

  Index AddNode( cFNumber aDfct ) override;

  void ChangeArc( cIndex name ,
		  cIndex nSS = Inf<Index>() , cIndex nEN = Inf<Index>() )
   override;

  void DelArc( cIndex name ) override;

  inline bool IsDeletedArc( cIndex name ) override;

  Index AddArc( cIndex Start , cIndex End , cFNumber aU , cCNumber aC )
   override; 

/*--------------------------------------------------------------------------*/
/*------------------------ SPECIALIZED INTERFACE ---------------------------*/
/*--------------------------------------------------------------------------*/

   void ShortestPathTree( void );

/**< Solver of the Shortest Path Tree Problem from the current Origin.
   (specified in the constructor or by SetOrigin(), see below)

   If LABEL_SETTING == 0, or if no Destination is speficied (Dst ==
   Inf<Index>() in SetDest() [see below]), the whole Shortest Path Tree (at
   least, the SPT of the component of the graph connected with Origin) is
   computed, otherwise the code stops as soon as the shortest path between
   Origin and Dest is computed.

   Note that methods such as MCFGetX(), MCFGetRC() and MCFGetFO() may need
   some complicate calculations in order to put the solution of the Shortest
   Path in the correct format; since these calculations change some of the
   internal data structures, it is not permitted to call again
   ShortestPathTree() after that any of these methods have been called. */

/*--------------------------------------------------------------------------*/

   inline void SetOrigin( cIndex NewOrg );

/**< Changes the Origin from which Shortest Paths are computed. */

/*--------------------------------------------------------------------------*/

   inline void SetDest( cIndex NewDst );

/**< Changes the Destination node of Shotest Paths. If LABEL_SETTING == 0, it
   has no influence since label correcting methods cannot stop before the
   whole SPT has been computed. Conversely, label setting algorithms can solve
   Origin-Dest Shortest Path Problems; therefore, it is possible to obtain
   shortest paths between Origin and a subset of the nodes, by calling
   ShortestPathTree() with one of the destinations, and controlling upon
   completion that all the desidered nodes have been visited (see Reached()
   below). If this is not the case, ShortestPathTree() can be invoked again
   with one of the unreached nodes, until they are all visited.

   If no Dest is given, or if Dest is set to Inf<Index>(), the whole Shortest
   Path Tree (at least, the SPT of the component of the graph connected with
   Origin) is computed. */

/*--------------------------------------------------------------------------*/

   void MCFGetX( Index ND , cIndex_Set DB , FRow F , Index_Set nms = 0 ,
		 cIndex strt = 0 , Index stp = Inf<Index>() );

/**< Like SPTree::MCFGetX( FRow , Index_Set , cIndex , Index ), except that
   the primal solution that is returned is relative only to the subset of
   destinations whose names are contained in the first ND entries of the
   vector DB.

   Note: node names in ND must be in 1 ... n irrespective of USENAME0. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   SPTree::FONumber MCFGetFO( Index ND , cIndex_Set DB );

/**< Like SPTree::MCFGetFO( void ), except that the cost that is returned is
   that of the primal solution relative only to the subset of destinations
   whose names are contained in the first ND entries of the vector DB.

   Note: node names in ND must be in 1 ... n irrespective of USENAME0. */

/*--------------------------------------------------------------------------*/

   inline bool Reached( cIndex i );

/**< Return true if a shortest path from Origin to i have already been
   computed; this can be used when LABEL_SETTING == 1 to determine if a
   shortest from Origin to i have been obtained as a by-product of the
   calculation of the shortest path between Origin and some other Dest. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set Predecessors( void );

/**< Return a cIndex* vector p[] such that p[ i ] is the predecessor of node
   i in the shortest path tree. If a node i has no predecessor, i.e.,
   i == Origin, i does not belong to the connected component of the origin or
   the computation have been stopped before reaching i, then p[ i ] == 0.

   \note if the name "0" is used for nodes, (USENAME0 == 1) then node names
         are internally "translated" of +1 to avoid it being used - the
	 the names reported in this vector will follow the same rule.

   For this reason, the first entry of p (*p) is not significative. */

/*--------------------------------------------------------------------------*/

   cIndex_Set ArcPredecessors( void );

/**< Return a cIndex* vector a[] such that a[ i ] is the index of the arc
   ( p[ i ] , i ), being p[] the vector returned by the above method, and
   with the same structure. If p[ i ] == 0, then a[ i ] is not significative:
   for the Origin (that has p[ Origin ] == 0), however, it is guaranteed that
   a[ Origin ] == Inf<Index>(). */

/*--------------------------------------------------------------------------*/

   inline Index Orig( void );

/**< Return the root of the SPT problem. */

/*--------------------------------------------------------------------------*/

   inline Index DestN( void );

/**< Return the number of destination nodes in the SPT problem. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   inline cIndex_Set Dests( void );

/**< Return the DestN()-vector containig the names of destination nodes in
   the SPT problem; the names are in increasing order and INF-terminated. */

/*--------------------------------------------------------------------------*/

   inline Index LenFS( cIndex i );

/**< Return the size of the Forward Star of node i. */

/*--------------------------------------------------------------------------*/

   inline Index ReadFS( cIndex i , cIndex h );

/**< Return the h-th arc in FS( i ) for h = 0, ... , LenFS( i ) - 1. */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~SPTree();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 struct FSElmnt {               // one entry of the Forward Star
                  CNumber Cst;  // cost of the arc
                  Index   Nde;  // end node of the arc
	          };

 typedef FSElmnt *FrwdStr;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

 Index Origin;       // the source
 Index Dest;         // the sink
 Index NDsts;        // total number of destinations;
 Index_Set DstBse;   // array of indices of the destinations

 Index_Set NdePrd;   // NdePrd[ i ] = predecessor of i in the shortest path
                     // NdePrd[ Origin ] = 0

 Index_Set ArcPrd;   // ArcPrd[ i ] = index of arc ( NdePrd[ i ] , i )
                     // ArcPrd[ Origin ] = 0

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
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void Initialize( void );

/* Initialize the data structures for a "cold start". */

/*--------------------------------------------------------------------------*/

   inline void ScanFS( cIndex mi );

/* Scans the Forward Star of mi and puts in Q those nodes whose distance
   label can be decreased by using an arc emanating from mi. */

/*--------------------------------------------------------------------------*/

   inline Index ExtractQ( void );

/* Extracts an element (depending on the particular algoritm) from the set Q:
   if Q is empty, returns 0. */

/*--------------------------------------------------------------------------*/

#if( ( SPT_ALGRTM == 0 ) || ( SPT_ALGRTM == 3 ) )

   inline void InsertQ( cIndex j );

#else

   inline void InsertQ( cIndex j , cCNumber label );

#endif

/* Inserts the node with name j and label label somewhere in Q: the label is
   not needed for LQueue and Djkstra algorithms. */

/*--------------------------------------------------------------------------*/

   inline void CalcArcP( void );

/* Calculates the ArcPrd[] vector. */

/*--------------------------------------------------------------------------*/

#if( ! SPT_STRTN )

   inline Index Startn( cIndex What );

/* Extract the starting node of the arc that is in position What in FS[]
   using a binary search on StartFS[]. */

#endif

/*--------------------------------------------------------------------------*/

 inline void MemAlloc( void );

 inline void MemDeAlloc( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 CRow Pi;            // node Potentials
 FRow B;             // node deficits vector

 FONumber FO;        // Objective Function value

 bool ReadyArcP;     // if the "arc predecessor" data structure has already
                     // been updated after a "final" ShortestPathTree() call 
 FrwdStr FS;         // the Forward Star (itself)

 Index_Set Q;        // the set of scanned nodes: Q[ i ] = INF ==> i \notin Q
 #if( SPT_ALGRTM <= 3 )
                     // here, Q is an array pointer implementation of a list,
                     // and *Q is the head of the list (node names are >= 1)
 #else
  Index_Set H;       // here, Q[ i ] tells the position of node i in the
		     // vector implementing the heap, and H is that vector
 #endif

 Index cFS;          // cardinality of the FS (m if DirSPT, 2m otherwise)
 Index tail;         // the tail element of the list, or the first free
                     // position in the heap

 // static members- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 static Index_Set Stack;  // temporary for flow calculation and
                          // reoptimization
 static Index InstCntr;   // counter of active instances
 static Index maxnmax;    // max value of nmax among all the instances

 #if( SPT_STRTN )
  #if( SAME_GRPH_SPT )
   static Index_Set Startn;
  #else
   Index_Set Startn;
  #endif
 #endif

 #if( SAME_GRPH_SPT && ( ! DYNMC_MCF_SPT ) )
  static Index_Set StrtFS;

  static Index_Set Dict;
  static Index_Set DictM1;

  static bool DirSPT;
 #else
  Index_Set StrtFS;  // position in FS[] where FS[ i ] begins
  #if( DYNMC_MCF_SPT )
   Index_Set LenFS;  // how many arcs there are in FS[ i ]
  #endif

  Index_Set Dict;    // arc dictionary: for each position in FS[], tells
                     // which arc is that one
  Index_Set DictM1;  // inverse of Dict: for each arc, tells where it stands
                     // in FS[] - if the graph is undirected, the two
		     // consecutive entries 2 * i and 2 * i + 1 tells the
		     // two positions of arc i in FS[]
  bool DirSPT;       // true if the graph is directed
 #endif

/*--------------------------------------------------------------------------*/

 };  // end( class SPTree )

/* @} end( group( SPTREE_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::MCFSNde( cIndex i )
{
 if( DirSPT ) 
  #if( SPT_STRTN )
   return( Startn[ i ] );
  #else
   return( Startn( DictM1[ i ] ) - USENAME0 );
  #endif
 else
  return( FS[ DictM1[ 2 * i + 1 ] ].Nde );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::MCFENde( cIndex i )
{
 if( DirSPT ) 
  return( FS[ DictM1[ i ] ].Nde - USENAME0 );
 else
  return( FS[ DictM1[ 2 * i ] ].Nde - USENAME0 );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber SPTree::MCFCost( cIndex i )
{
 if( DirSPT ) 
  return( FS[ DictM1[ i ] ].Cst );
 else
  return( FS[ DictM1[ 2 * i ] ].Cst );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::FNumber SPTree::MCFUCap( cIndex i )
{
 #if( ! DYNMC_MCF_SPT )
  cIndex pos = DictM1[ i ];
  #if( SAME_GRPH_SPT )
   if( FS[ pos ].Cst == Inf<CNumber>() )
  #else
   if( pos >= StrtFS[ n + 1 ] )
  #endif
    return( 0 );
   else
 #endif
    return( - B[ Origin ] );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::FNumber SPTree::MCFDfct( cIndex i )
{
 return( B[ i + 1 ] );
 }

/*--------------------------------------------------------------------------*/

inline bool SPTree::IsClosedArc( cIndex name )
{
 #if( DYNMC_MCF_SPT )
  cIndex pos = DictM1[ name ];     // current position of arc name
  #if( SPT_STRTN )
   cIndex nde = Startn[ name ];    // start node of arc name
  #else
   cIndex nde = Startn( pos );     // start node of arc name
  #endif

  return( pos < StrtFS[ nde ] + LenFS[ nde ] );
 #else
  return( false );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool SPTree::IsDeletedArc( cIndex name )
{
 return( SPTree::IsClosedArc( name ) );
 }

/*--------------------------------------------------------------------------*/

inline void SPTree::SetOrigin( cIndex NewOrg )
{
 if( Origin != NewOrg + USENAME0 ) {
  Origin = NewOrg + USENAME0;
  status = MCFClass::kUnSolved;
  }
 }

/*--------------------------------------------------------------------------*/

inline void SPTree::SetDest( cIndex NewDst )
{
 if( Dest != NewDst + USENAME0 ) {
  #if( LABEL_SETTING )
   Dest = NewDst + USENAME0;
  #endif
  status = MCFClass::kUnSolved;
  }
 }

/*--------------------------------------------------------------------------*/

inline bool SPTree::Reached( cIndex i )
{
 return( ( Pi[ i ] < Inf<CNumber>() ) && ( Q[ i ] == Inf<Index>() ) );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::cIndex_Set SPTree::Predecessors( void )
{
 return( NdePrd );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::Orig( void )
{
 return( Origin );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::DestN( void )
{
 return( NDsts );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::cIndex_Set SPTree::Dests( void )
{
 return( DstBse );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::LenFS( cIndex i )
{
 #if( DYNMC_MCF_SPT )
  return( LenFS[ i ] );
 #else
  return( StrtFS[ i + 1 ] - StrtFS[ i ] );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::Index SPTree::ReadFS( cIndex i , cIndex h )
{
 return( Dict[ StrtFS[ i ] + h ] );
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SPTree.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File SPTree.h -----------------------------*/
/*--------------------------------------------------------------------------*/
