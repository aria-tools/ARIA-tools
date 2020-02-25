/*--------------------------------------------------------------------------*/
/*------------------------- File RelaxIV.h ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Linear Min Cost Flow problems solver, based on the RELAXIV code by
 * D. Bertsekas and P. Tseng, as described in
 *
 *    Bertsekas, Dimitri P., and Paul Tseng.
 *    "RELAX-IV: A faster version of the RELAX code for solving minimum
 *     cost flow problems." (1994), Report LIDS-P-2276, MIT.
 *
 * Conforms to the standard (MCF) interface defined in MCFClass.h.
 *
 * RelaxIV is based on a primal-dual algorithm which essentially operates
 * as follows: a pseudoflow (a flow vector which satisfies bound and
 * non-negativity constraints but not necessarily flow conservation
 * constraints) is kept which satisfies complementarity slackness conditions
 * with the current vector of potentials; that is, only the flow on arcs
 * whose reduced cost
 * \f[
 *  RC[ i , j ] =  C[ i , j ] - Pi[ j ] + Pi[ i ]
 * \f]
 * is zero can be chosen to any value between 0 and the capacity, while arcs
 * with reduced cost < 0 are saturated (fixed to their capacity) and arcs
 * with reduced cost > 0 are empty (fixed to 0).
 *
 * The algorithm attempts to convert the pseudoflow into a flow (i.e.,
 * to satisfy the flow conservation constraints) by essentially running
 * a max-flow algorithm of the augmenting path type. If the flow is found then
 * this is an optimal solution of the problem and the algorithm is stopped.
 * Otherwise, a saturated cut is identified which separates the origins (nodes
 * not yet producing enough flow) to the destinations (nodes not yet consuming
 * enough flow); this cut is used to modify the potentials, thereby creating
 * new arcs with zero reduced cost, which can be used to push further flow
 * from the origins to the destinations. If no such arcs can be created the
 * problem is declared unfeasible. Much care is devoted to stop the max-flow
 * computation as soon as a proof that the set of potentials is not optimal,
 * in order to reach as soon as possible a dual optimal solution, and to
 * re-use all available information to "warm start" the max-flow computation
 * after a change in the potentials.
 *
 * \warning The original code has been written for integer data only.
 *          By properly setting the flow and cost tolerances [see
 *          SetEps****() in MCFClass.h] we have always been able to solve
 *          any MCF that we could throw at the solver, but in principle this
 *          kind of algorithm may fail to converge with nonintegral data, so
 *          consider yourselves warned.
 *
 * \warning A private type SIndex is defined which is intended to hold arc
 *          and node indices "with a sign", used to represent orientation.
 *          This has to be "in sync" with Index, in the sense that for every
 *          unsigned index value in Index, the two signed values should be
 *          feasible in SIndex. In other words, either Index is not using at
 *          least half of its feasible values, or SIndex has to be a "bigger"
 *          data type than Index. The default value for SIndex is int.
 *
 * \version 1.84
 *
 * \date 16 - 09 - 2019
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Dimitri P. Bertsekas \n
 *         Lab. for Information and Decision Systems \n
 *         Massachusetts Institute of Technology \n
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Paul Tseng \n
 *         Department of Mathematics \n
 *         University of Washington \m
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * Copyright &copy 1996 - 2019 by Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __RelaxIV
 #define __RelaxIV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup RELAXIV_MACROS Compile-time switches in RelaxIV.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

/*----------------------------- SAME_GRPH_RIV ------------------------------*/

#define SAME_GRPH_RIV 0

/**< Decides if all MCFClass instances share the same graph.
   If SAME_GRPH_RIV > 0, then all the instances of the class will work on the
   same "topological" network, while the costs, capacities and supplies can
   change from one instance to another. This allows implementations to share
   some data structures describing the graph, e.g. by declaring them "static",
   saving memory when multiple instances of the solver are active at the same
   time. However, this also obviously inhibits some type of changes in the
   topology of the graph [see DYNMC_MCF_RIV below]. */

/*----------------------------- DYNMC_MCF_RIV ------------------------------*/

#define DYNMC_MCF_RIV 3

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_RIV > 0, the methods of the public interface of class that
   allow to change the topology of the underlying network are actually
   implemented. Possible values of this macro are:

   - 0 => the topology of the graph cannot be changed;

   - 1 => the methods that "close" arcs and delete nodes are implemented;

   - 2 => the methods that "open" previously closed arcs and add nodes are
          implemented;

   - 3 => the methods that change the start and end node of a (possibly
          "closed") arc, delete and create new arcs are implemented.

   As long as DYNMC_MCF_RIV < 3, SAME_GRPH_RIV can possibly be > 0; less and
   less data structures can be shared as DYNMC_MCF_RIV increases.
   DYNMC_MCF_RIV == 3 implies SAME_GRPH_RIV == 0. */

/*-------------------------------- AUCTION ---------------------------------*/

#define AUCTION 0

/**< Decides if the auction/shortest paths inizialization procedure is used.
   If AUCTION == 1, then an auction/shortest paths inizialization procedure
   is provided [see SetPar() below] that has been reported to make the
   RelaxIV algorithm run faster on some classes of instances.
   The auction initialization essentially "spreads" the imbalances around
   the graph by performing some steps of the "pure" epsilon-relaxation
   method: this should produce "short" augmenting steps, that seem to be
   the best situation for RelaxIV.

   By setting AUCTION == 0, some memory is saved. */

/*-------------------------- RELAXIV_STATISTICS ----------------------------*/

#define RELAXIV_STATISTICS 0

/**< If RELAXIV_STATISTICS > 0, then statistic information about the behaviour
   of the Relaxation algorithm is computed. */

/*@} -----------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS RelaxIV --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup RELAXIV_CLASSES Classes in RelaxIV.h
    @{ */

/** The RelaxIV class derives from the abstract base class MCFClass, thus
    sharing its (standard) interface, and implements a Relaxation algorithm
    for solving (Linear) Min Cost Flow problems. */

class RelaxIV : public MCFClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Public enum describing the possible parameters of the MCF solver,
    "extended" from MCFClass::MCFParam, to be used with the methods
    SetPar() and GetPar(). */

  enum MCFRParam { kAuction = kLastParam     ///< crash initialization
                   };


/*--------------------------------------------------------------------------*/
/** Public enum describing the more file formats in RelaxIV::WriteMCF(). */

  enum RIVFlFrmt { kCLP = kMPS + 1 ,  ///< the "LP" format
		   kRIV               ///< RelaxIV-specific format
                   };

/*--------------------------------------------------------------------------*/

  typedef bool           *Bool_Vec;        ///< vector of booleans

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   RelaxIV( cIndex nmx = 0 , cIndex mmx = 0 );

/**< Constructor of the class, as in MCFClass::MCFClass(). */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
		 cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
		 cIndex_Set pEn = NULL ) override;

/**< Inputs a new network, as in MCFClass::LoadNet().

   Arcs with pC[ i ] == Inf<CNumber>() do not "exist". If DYNMC_MCF_RIV > 0,
   these arcs are "closed".

   If DYNMC_MCF_RIV == 0 but SAME_GRPH_RIV > 0, these arcs are dealt with
   exactly like as if they had pU[ i ] == 0, i.e., as "normal" arcs with
   zero capacity. These arcs can be put back into the formulation by simply
   changing their capacity (and cost). Note that, however, this is less
   efficient than eliminating them explicitly from the problem.

   If DYNMC_MCF_RIV == 0 and SAME_GRPH_RIV == 0, these arcs are just removed
   from the formulation. However, they have some sort of a "special status"
   (after all, if the user wants to remove them completely he/she can just
   change the data), in that they are still counted into the number of arcs
   of the graph and they will always have 0 flow and Inf<CNumber>()
   reduced cost as "closed" or "deleted" arcs.

   If SAME_GRPH_RIV == 1 (==> DYNMC_MCF_RIV < 3), pSn and pEn passed to any
   instance after the first that is constructed are ignored, thus,
   pSn == pEn == NULL is allowed in this case. Note that pn and pm are
   *not* ignored. */

/*--------------------------------------------------------------------------*/
/// set integer parameters of the algorithm
/** Set integer parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - kAuction: if set to kYes, the auction/shortest paths initialization is
               used in SolveMCF() to generate the starting solution; if
	       set to kNo (default), then the default initialization based on
	       special single-node relaxation iterations is used instead.
	       Note that this parameter is *ignored* if AUCTION == 0. */

   virtual void SetPar( int par , int val ) override
   {
    if( par == kAuction ) {
     #if( AUCTION )
      crash = ( val == kYes ) ? TRUE : FALSE;
     #else
      if( val == kYes )
       throw( MCFException( "Auction initialization not available" ) );
     #endif
     }
    else
     MCFClass::SetPar( par , val );
  }

/*--------------------------------------------------------------------------*/
/// set double parameters of the algorithm
/** Set double parameters of the algorithm. This only calls the base class
 * method. It should not be necessary, but sometimes it is. */

   virtual void SetPar( int par , double val ) override
   {
    MCFClass::SetPar( par , val );
    }

/*--------------------------------------------------------------------------*/

   virtual inline void GetPar( int par , int &val ) const override
   {
    if( par == kAuction )
     #if( AUCTION )
      val = crash ? kYes : kNo;
     #else
      val = kNo;
     #endif
    else
     MCFClass::GetPar( par , val );
    }

/**< This method returns one of the integer parameter of the algorithm.

   @param par  is the parameter to return [see SetPar( int ) for comments];

   @param val  upon return, it will contain the value of the parameter.

   Apart from the parameters of the base class, this method handles kAuction.
   */

/*--------------------------------------------------------------------------*/

 virtual inline void GetPar( int par , double &val ) const override
 {
  MCFClass::GetPar( par , val );
  }

 /*--------------------------------------------------------------------------*/

   void PreProcess( void ) override;

/**< If this method is called, a preprocessing phase is performed trying to
   reduce the arc capacities. This may sometimes help in speeding up the
   solution of the problem, but may also change the capacities returned by
   MCFUCap[s]() [see below].

   For this method to work properly, arc capacities, node deficits and the
   topology of the graph must have already been provided with LoadNet()
   [see above].

   This method can be called more than once, for instance whenever the
   capacities of some arcs or the deficits of some nodes are changed;
   however, it destroys the provious optimal solution (if any), forcing the
   algorithm to restart from scratch. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   void MCFGetX( FRow F , Index_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   cFRow MCFGetX( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFGetRC( CRow CR , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   cCRow MCFGetRC( void ) override;

   CNumber MCFGetRC( cIndex i ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFGetPi( CRow P , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

/**< Same meaning as MCFClass::MCFGetPi().

   \note if both AUCTION == 0 and DYNMC_MCF_RIV <= 1, no internal memory for
   the vector of potentials is allocated; hence, even if nms != NULL
   MCFGetPi( P ) first constructs the full vector of potentials in P and
   then selects only the components in nms. Therefore, memory can be written
   even after the | nms |-th element of P. For the same reason, memory can be
   written even after themore the (stp - strt)-th element of P. */

   cCRow MCFGetPi( void ) override;

/**< Same meaning as MCFClass::MCFGetPi().

   \note if MCFGetPi( void ) returns a pointer, this is a pointer to a static
   shared data structure that can be corrupted by calls to methods (e.g.
   MCFGetPi() itself) for other active instances of RelaxIV. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   RelaxIV::FONumber MCFGetFO( void ) override;

/*--------------------------------------------------------------------------*/

   MCFClass::MCFStatePtr MCFGetState( void ) override;

/**< Same meaning as MCFClass::MCFGetState().

   The state of the algorithm is the pair S = ( X[] , RC[] ) of the arc
   flows and reduced costs. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFPutState( MCFClass::MCFStatePtr S ) override;

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline Index MCFSNde( cIndex i ) override;

   inline Index MCFENde( cIndex i ) override;

   inline cIndex_Set MCFSNdes( void ) override;

/**< Same meaning as MCFClass::MCFSNdes().

   \note MCFSNdes() returns a pointers to a (read-only) vector containing
         the arc start nodes *only if USENAME0 == 0*; otherwise, it returns
	 NULL. */

   inline cIndex_Set MCFENdes( void ) override;

/**< Same meaning as MCFClass::MCFENdes().

   \note MCFENdes() returns a pointers to a (read-only) vector containing
         the arc end nodes *only if USENAME0 == 0*; otherwise, it returns
	 NULL. */

/*--------------------------------------------------------------------------*/

   void MCFCosts( CRow Costv , cIndex_Set nms = NULL  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline CNumber MCFCost( cIndex i ) override;

   inline cCRow MCFCosts( void ) override;

/*--------------------------------------------------------------------------*/

   void MCFUCaps( FRow UCapv , cIndex_Set nms = NULL  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline FNumber MCFUCap( cIndex i ) override;

   inline cFRow MCFUCaps( void ) override;

/*--------------------------------------------------------------------------*/

   void MCFDfcts( FRow Dfctv , cIndex_Set nms = NULL  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   inline FNumber MCFDfct( cIndex i ) override;

   inline cFRow MCFDfcts( void ) override;

/*--------------------------------------------------------------------------*/

   void WriteMCF( ostream &oStrm , int frmt = 0 ) override;

/**< Extends MCFClass::WriteMCF() to support two new formats:

   - kCLP is the "LP" format read by several LP solvers;

   - kRIV is the following RelaxIV-specific format:

          - < number of nodes > < number of arcs >

	  - for( < each arc > )
	    < start node > < end node > < reduced_capacity > < reduced_cost >

	  - for( < each node > )
	    < reduced flow deficit at node >

	  \note the data of the problem in this format is not that of the
	        original problem, but rather that of the "reduced" problem
	        corresponding to the current pair (flow, potential) of the
		relaxation algorithm. */

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgCost( Index arc , cCNumber NCost ) override;

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgDfct( Index nod , cFNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) override;

   void ChgUCap( Index arc , cFNumber NCap  ) override;

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

   inline int MCFiter();
   ///< Total number of (single-node or multinode) iterations

   inline int MCFaug();
   ///< Number of flow augmentations

 #if( RELAXIV_STATISTICS )
   inline int MCFmulti();
   ///< Number of multinode iterations

   inline int MCFascnt();
   ///< Number of dual ascent steps

  #if( AUCTION )
   inline int MCFauct();  
   ///< Number of iterations in the Auction() initialization
  #endif
 #endif

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  virtual ~RelaxIV();

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

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
/*---------------------------- PRIVATE TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

  typedef int		  SIndex;            ///< an index with a sign
  typedef SIndex          *SIndex_Set;       ///< set (array) of SIndex
  typedef const SIndex    cSIndex;           ///< a read-only SIndex
  typedef cSIndex        *cSIndex_Set;       ///< read-only SIndex array

/*--------------------------------------------------------------------------*/

   class RIVState : public MCFClass::MCFState {
    public:

     RIVState( cIndex m );
     ~RIVState();

     FRow Flow;
     CRow RedCost;
     };

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- called in SolveMCF() ---------------------------*/
/*--------------------------------------------------------------------------*/

 inline void init_tree();

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void init_standard( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline FNumber svblncdarcs( cIndex node ,
			     cIndex_Set tfst1 , cIndex_Set tnxt1 ,
			     cIndex_Set tfst2 , cIndex_Set tnxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline FNumber dascnt( cIndex node , CNumber &delprc , cIndex_Set F1 ,
		        cIndex_Set Nxt1 , cIndex_Set F2 , cIndex_Set Nxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void relist( cIndex node );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void AugFlow( cIndex augnod , cIndex root , cIndex node_p ,
		      cIndex node_n , cIndex_Set Term1 , cIndex_Set Term2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline bool Ascnt( cFNumber sdm , FNumber delx , Index &nlabel ,
		    bool &Switch , Index &nscan , Index &curnode ,
		    cIndex_Set Term1 , cIndex_Set Term2 ,  cIndex_Set F1 ,
		    cIndex_Set Nxt1 , cIndex_Set F2 , cIndex_Set Nxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 #if( AUCTION )
  inline void auction( void );
 #endif

/*--------------------------------------------------------------------------*/
/*----------------------- called in init_standard --------------------------*/
/*--------------------------------------------------------------------------*/

 inline CNumber nxtbrkpt( cIndex_Set t_St1 , cIndex_Set NSt1 ,
			  cIndex_Set t_St2 , cIndex_Set NSt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline CNumber mvflw1( cIndex arc , FRow tDfct , FRow tDDNeg ,
		        cIndex_Set Term , FRow Flow1 , FRow Flow2 );

 inline CNumber mvflw2( cIndex arc , FRow tDfct , FRow tDDPos ,
		        cIndex_Set Term , FRow Flow1 , FRow Flow2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void decrsRC( cIndex arc , CNumber trc , cCNumber delprc ,
		      CNumber &nxtbrk , FRow tDD1 , FRow DD2 ,
		      cIndex_Set Term );

 inline void incrsRC( cIndex arc , CNumber trc , cCNumber delprc ,
		      CNumber& nxtbrk , FRow tDD1 , FRow DD2 ,
		      cIndex_Set Term );

/*--------------------------------------------------------------------------*/
/*--------------------------- called in Chg**** ----------------------------*/
/*--------------------------------------------------------------------------*/

 inline void chgcsti( cIndex i , CNumber NCost );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void chgcapi( cIndex i , cFNumber NCap );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( DYNMC_MCF_RIV )

 inline void delarci( cIndex arc );

 #if( DYNMC_MCF_RIV > 1 )

  inline void addarci( cIndex arc );

 #endif
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline void cmptprices( void );

/*--------------------------------------------------------------------------*/

 inline void MemAlloc( void );

 inline void MemDeAlloc( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 FRow X;                  // arc Flows
 FRow U;                  // arc residual capacities
 FRow Cap;                // arc Capacities

 CRow C;                  // arc Costs
 CRow RC;                 // arc Reduced Costs

 FRow B;                  // node deficits vector
 FRow Dfct;               // node residual deficits

 FONumber FO;             // Objective Function value

 Index_Set tfstou;        // tfstou e tnxtou describe the subsets of the
 Index_Set tfstin;        // forward stars composed by balanced arcs;
 Index_Set tnxtin;        // tfstin e tnxtin have the same function but for
 Index_Set tnxtou;        // backward stars

 Index nb_pos;            // number of "directed" ...
 Index nb_neg;            // ... and "inverse" balanced arcs

 #if( DYNMC_MCF_RIV > 2 )
  Index ffp;              // first free position in arc vectors
 #endif

 #if( AUCTION )
  bool crash;             // true => initialization is perfomed by the
                          // auction routine, false => it is performed by
			  // single node relaxation iterations
 #endif

 int iter;                // number of iterations (of both types)
 int num_augm;            // number of flow augmentation steps
 #if( RELAXIV_STATISTICS )
  int nmultinode;         // number of multinode iterations
  int num_ascnt;          // number of multinode ascent steps
  #if( AUCTION )
   int nsp;               // n. of auction/shortest path iterations
  #endif
 #endif

 Index error_node;        // error handling variables: if the problem is
 Index error_info;        // found to be unfeasible, these variables contain
                          // a description of the kind of unfeasibility and
                          // where it is located. error_info is
                          // 1 unfeasibility detected in PreProcessing(): out
                          //   capacity of node error_node < - deficit
                          // 2 unfeasibility detected in PreProcessing(): in
                          //   capacity of node error_node < deficit
                          // 3 exit during initialization by single node
                          //   iterations: dual ascent feasible ray was found
                          //   while increasing price of node error_node;
                          // 4 exit during initialization by single node
                          //   iterations: dual ascent feasible ray was found
                          //   while decreasing price of node error_node;
                          // 5 dual ascent feasible ray was found during a
                          //   relaxation iterazion at node error_node with
                          //   positive deficit;
                          // 6 dual ascent feasible ray was found during a
                          //   relaxation iterazion at node error_node with
                          //   negative deficit;
                          // 7 dual ascent feasible ray was found during a
                          //   multinode relaxation iteration, error_node is
                          //   the starting node of the iteration;
                          // 8 problem has been detected unfeasible in
                          //   Auction() initialization.

 // potentially static members, depending on SAME_GRPH_RIV - - - - - - - - - -

 #if( SAME_GRPH_RIV )
  static Index_Set Startn;  // Start ...
  static Index_Set Endn;    // .. and End node of each edge
 #else
  Index_Set Startn;
  Index_Set Endn;
 #endif

 #if( SAME_GRPH_RIV && ( ! DYNMC_MCF_RIV ) )
  static Index_Set FOu;     // index of the first edge exiting from node
  static Index_Set FIn;     // index of the first edge entering into node
  static Index_Set NxtOu;   // for each edge a, NxtOu[ a ] is the next edge 
			    // exiting from Startn[ a ]
  static Index_Set NxtIn;   // analogous for entering arcs
 #else
  Index_Set FOu;
  Index_Set FIn;
  Index_Set NxtOu;
  Index_Set NxtIn;
 #endif

 #if( AUCTION || ( DYNMC_MCF_RIV > 1 ) )
  static RelaxIV *PiOwnr;  // the instance who calculated Pi the last time
                           // (NULL if they have not been computed)
 #endif

 // static members- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 static CRow Pi;          // node Potentials

 static Bool_Vec mark;    // various temporaries for multinode iterations
 static Index_Set save;
 static Index_Set label;
 static SIndex_Set Prdcsr;

 static Bool_Vec scan;    // in multinode iteration, scan denote that a
		          // node belongs to S in the cut
 static Index_Set queue;  // queue of non zero deficit nodes
 static Index lastq;      // index of the last element in the queue
 static Index prvnde;     // index of the element preceding lastqueue

 static FRow DDNeg;       // positive directional derivative at nodes
 static FRow DDPos;       // negative directional derivative at nodes

 #if( AUCTION )
  static CRow SB_level;   // various temporaries used in Auction()
  static SIndex_Set extend_arc;
  static SIndex_Set SB_arc;
  static Index_Set FpushF;
  static Index_Set NxtpushF;
  static Index_Set FpushB;
  static Index_Set NxtpushB;
 #endif

 static Index InstCntr;   // counter of active instances
 static Index maxnmax;    // max value of nmax among all the instances
 static Index maxmmax;    // max value of nmax among all the instances

/*--------------------------------------------------------------------------*/

 };  // end( class RelaxIV )

/* @} end( group( RELAXIV_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline MCFClass::cIndex_Set RelaxIV::MCFSNdes( void )
{
 #if( USENAME0 )
  return( NULL );
 #else
  return( Startn + 1 );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::cIndex_Set RelaxIV::MCFENdes( void )
{
 #if( USENAME0 )
  return( NULL );
 #else
  return( Endn + 1 );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::Index RelaxIV::MCFSNde( cIndex i )
{
 return( Startn[ i + 1 ] - USENAME0 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::Index RelaxIV::MCFENde( cIndex i )
{
 return( Endn[ i + 1 ] - USENAME0 );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::cCRow RelaxIV::MCFCosts( void )
{
 return( C + 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::CNumber RelaxIV::MCFCost( cIndex i )
{
 return( C[ i + 1 ] );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::cFRow RelaxIV::MCFUCaps( void )
{
 return( Cap + 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::FNumber RelaxIV::MCFUCap( cIndex i )
{
 return( Cap[ i + 1 ] );
 }

/*--------------------------------------------------------------------------*/

inline MCFClass::cFRow RelaxIV::MCFDfcts( void )
{
 return( B + 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

inline MCFClass::FNumber RelaxIV::MCFDfct( cIndex i )
{
 return( B[ i + 1 ] );
 }

/*--------------------------------------------------------------------------*/

inline bool RelaxIV::IsClosedArc( cIndex name )
{
 #if( DYNMC_MCF_RIV > 2 )
  return( ( RC[ name + 1 ] == Inf<CNumber>() ) &&
	  ( Startn[ name + 1 ] < Inf<Index>() ) );
 #elif( DYNMC_MCF_RIV )
  return( RC[ name + 1 ] == Inf<CNumber>() );
 #else
  return( false );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool RelaxIV::IsDeletedArc( cIndex name )
{
 #if( DYNMC_MCF_RIV > 2 )
  return( Startn[ name + 1 ] == Inf<Index>() );
 #else
  return( false );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline int RelaxIV::MCFiter( void )
{
 return( iter );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline int RelaxIV::MCFaug( void ) 
{
 return( num_augm );
 }

/*--------------------------------------------------------------------------*/

#if( RELAXIV_STATISTICS )
 inline int RelaxIV::MCFmulti( void )
 {
  return( nmultinode );
  }
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 inline int RelaxIV::MCFascnt( void )
 {
  return( num_ascnt );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 #if( AUCTION )
  inline int RelaxIV::MCFauct( void )
  {
   return( nsp );
   }
 #endif
#endif

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* RelaxIV.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File RelaxIV.h ----------------------------*/
/*--------------------------------------------------------------------------*/
