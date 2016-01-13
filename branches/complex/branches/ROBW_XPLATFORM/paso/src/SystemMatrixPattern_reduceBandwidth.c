/*  Finley_Reduce determines a row and column permutation which, when applied to a given sparse matrix, produces a permuted */
/*  matrix with a smaller bandwidth and profile. the input array is a connection table which represents the                 */
/*  indices of the nonzero elements of the matrix, a.  the algorithm is described in terms of the adjacency graph which     */
/*  has the characteristic that there is an edge (connection) between nodes i and j if a(i,j) .ne. 0 and i .ne. j.          */


/*    NDSTK(NR,D1)        D1 IS >= MAXIMUM DEGREE OF ALL NODES.
/*    Label(D2)            D2 AND NR ARE >= THE TOTAL NUMBER OF
/*    NewLabel(D2+1)         NODES IN THE GRAPH.
/*    Degree(D2)            STORAGE REQUIREMENTS CAN BE SIGNifICANTLY

/*    AssignedLevel(D2)             DECREASED FOR IBM 360 AND 370 COMPUTERS
/*    LevelToNode(D2)           BY REPLACING INTEGER NDSTK BY
/*    NodeToLevel(D2)           INTEGER*2 NDSTK IN SUBROUTINES Finley_Reduce,
/*    ConnectedComponents(D2)          Finley_Reduce_setDegree, Finley_Reduce_findDiameter, Finley_Reduce_dropTree AND NUMBER.
/*  COMMON INFORMATION--THE FOLLOWING COMMON BLOCK MUST BE IN THE
/*  CALLING ROUTINE.
/*    COMMON/GRA/N,IDPTH,ndeg_max
/*  EXPLANATION OF INPUT VARIABLES--
/*    NDSTK-     CONNECTION TABLE REPRESENTING GRAPH.
/*               pattern->index[iptr]=NODE NUMBER OF JTH CONNECTION TO NODE
/*               NUMBER I.  A CONNECTION OF A NODE TO ITSELF IS NOT
/*               LISTED.  EXTRA POSITIONS MUST HAVE ZERO FILL.
/*    NR-        ROW DIMENSION ASSIGNED NDSTK IN CALLING PROGRAM.
/*    Label[i]-   NUMBERING OF ITH NODE UPON INPUT. if NO NUMBERING EXISTS THEN Label[i]=I.
/*    N-         NUMBER OF NODES IN GRAPH (EQUAL TO ORDER OF MATRIX).
/*    ndeg_max-      MAXIMUM DEGREE OF ANY NODE IN THE GRAPH.
/*    Degree[i]-   THE DEGREE OF THE ITH NODE.
/*  EXPLANATION OF OUTPUT VARIABLES--
/*    NewLabel[i]-  THE NEW NUMBER FOR THE ITH NODE.
/*    newBandwidth-      THE BANDWIDTH AFTER numbering .
/*    newProfile-      THE PROFILE AFTER numbering .
/*    IDPTH-     NUMBER OF LEVELS IN Finley_Reduce LEVEL STRUCTURE.
/*                                                                                                       */
/*  The following only have meaning if the graph was connected:

/*    AssignedLevel[i]-    INDEX INTO LevelToNode TO THE FIRST NODE IN LEVEL I. AssignedLevel(I+1)-AssignedLevel[i]= number of nodes in ith leveL
/*    LevelToNode-     node numbers listed by level
/*    NodeToLevel[i]-  the level assigned to node i by Finley_Reduce.

/*  WORKING STORAGE VARIABLE--
/*    ConnectedComponents
/*  LOCAL STORAGE--
/*    COMMON/CC/-SUBROUTINES Finley_Reduce, SORT2 AND PIKAssignedLevel ASSUME THAT
/*               THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
/*               SUBROUTINE Finley_Reduce_findDiameter ASSUMES THAT THERE ARE AT MOST
/*               100 NODES IN THE LAST LEVEL.
/*    COMMON/AssignedLevelW/-SUBROUTINES Finley_Reduce_setup AND PIKAssignedLevel ASSUME THAT THERE
/*               ARE AT MOST 100 LEVELS.
      SUBROUTINE Finley_Reduce(NDSTK, NR, Label, NewLabel, Degree, AssignedLevel, LevelToNode,NodeToLevel, ConnectedComponents, newBandwidth, newProfile, NIN, ndeg_maxIN)
/* USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
	integer NIN, ndeg_maxIN
      INTEGER NDSTK
      INTEGER StartNode, ReverseNode, NewLabel, ConnectedComponentCounter, SORT2, LastAvailableNodeNumber, ConnectedComponents,
     * SIZE, STPT, FirstAvailableNodeNumber
      COMMON /GRA/ N, IDPTH, ndeg_max
/* IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
      COMMON /CC/ ConnectedComponentCounter, SIZE(50), STPT(50)
      COMMON /AssignedLevelW/ NHIGH(100), NLOW(100), NACUM(100)
      DIMENSION ConnectedComponents(1), Label(1)
      DIMENSION NDSTK(NR,1), AssignedLevel(1), LevelToNode(1), NodeToLevel(1), NewLabel(1),
     * Degree(1)
      newBandwidth = 0
      newProfile = 0
/* SET NewLabel[i]=0 FOR ALL I TO INDICATE NODE I IS UNNUMBERED
c
c	N = NR
c	do i=0,i<N,++i)
c		if(ndeg_max<Degree[i]) then
c		ndeg_max = Degree[i]
c		end if
c	end do
c

 set Degree 
     N=NIN
     ndeg_max = ndeg_maxIN
     for (i=0,i<N,++i) {
        Degree[i]=ptr[i+1]-ptr[i];
        NewLabel[i] = 0;
     }

    /* COMPUTE DEGREE OF EACH NODE AND ORIGINAL BANDWIDTH AND PROFILE*/

    CALL Finley_Reduce_setDegree(NDSTK, NR, Degree, Label, initialBandwidth, initialProfile)

    /* FirstAvailableNodeNumber = low end of available numbers for numbering  */
    /* LastAvailableNodeNumber = high end of available numbers for numbering  */
    FirstAvailableNodeNumber = 1;
    LastAvailableNodeNumber = N;

    /* number the nodes of degree zero */

    for (i=0,i<N,++i) {
        if (Degree[i]==0) {
           NewLabel[i] = LastAvailableNodeNumber;
           LastAvailableNodeNumber--;
        }
   }
/* find an unnumbered node of min degree to start on */
   while (FirstAvailableNodeNumber<=LastAvailableNodeNumber) {

      LowestDegree = ndeg_max + 1;
      Flag = 1;
      SDIR = 1;
      for (i=0,i<N,++i) {
        if (Degree[i]<LowestDegree && NewLabel[i]<=0) {
           LowestDegree = Degree[i];
           StartNode = i;
        }
      }

/* find pseudo-diameter and associated level structures: */

/* StartNode and ReverseNode are the ends of the diameter and LevelToNode and NodeToLevel are the respective level structures. */

      CALL Finley_Reduce_findDiameter(StartNode, ReverseNode, NDSTK, NR, Degree, AssignedLevel, LevelToNode,NodeToLevel, ConnectedComponents, IDFLT)
      /* Flag INDICATES THE END TO BEGIN NUMBERING ON */
      if (Degree(StartNode)>Degree(ReverseNode)) {
         Flag = -1;
         StartNode = ReverseNode;
      }
      CALL Finley_Reduce_setup(AssignedLevel, LevelToNode, NodeToLevel)
      /* find all the connected components  (ConnectedComponentCounter counts them) */
      ConnectedComponentCounter = 0;
      LROOT = 1
      available_NodesInTree = 1
      for (i=0,i<N,++i) {
        if (AssignedLevel[i]==0) {
           ConnectedComponentCounter++;
           STPT(ConnectedComponentCounter) = LROOT
           CALL Finley_Reduce_dropTree(I, NDSTK, NR, AssignedLevel, ConnectedComponents, Degree, LastLevelWidth, BottomLevel,available_NodesInTree, max_LevelWidth, N)
           SIZE(ConnectedComponentCounter) = BottomLevel + LastLevelWidth - LROOT
           LROOT = BottomLevel + LastLevelWidth
           available_NodesInTree = LROOT
        }
      }
/* on return from PIKAssignedLevel, DirectionLargestComponent indicates the direction the largest component fell. */
/* DirectionLargestComponent is modIFied now to indicate the numbering direction. num is set to the proper value for this direction. */
      if (SORT2(DMY)!=0) PIKAssignedLevel(LevelToNode, NodeToLevel, ConnectedComponents, IDFLT, DirectionLargestComponent)
      DirectionLargestComponent = DirectionLargestComponent*Flag;
      NUM= (DirectionLargestComponent<0) ?  LastAvailableNodeNumber :FirstAvailableNodeNumber;

      CALL NUMBER(StartNode, NUM, NDSTK, NodeToLevel, Degree, NewLabel, LevelToNode,
     * AssignedLevel, NR, Flag, newBandwidth, newProfile, ConnectedComponents, DirectionLargestComponent)

      /* update LastAvailableNodeNumber or FirstAvailableNodeNumber after numbering */

      if (DirectionLargestComponent<0) LastAvailableNodeNumber = NUM;
      if (DirectionLargestComponent>0) FirstAvailableNodeNumber = NUM;

   }

   /* if original numbering is better than new one, set up to return it */

   if (newBandwidth > initialBandwidth) {
         for (i=0,i<N,++i)
             NewLabel[i] = Label[i];
         *newBandwidth = initialBandwidth;
         *newProfile = initialProfile;
   }
}

/*  Finley_Reduce_setDegree computes the degree of each node and stores it in the array Degree. */
/*  the bandwidth and profile for the original or input numbering  of the graph is computed also. */

void Finley_Reduce_setDegree(  pattern, dim_t *Degree,dim_t* Label,dim_t *bandwidth,dim_t *initialProfile) {

      *bandwidth = 0;
      *profile = 0;
      for (i=0,i<pattern->N,++i) {
        Degree[i] = 0;
        max_diff = 0;
        for (iptr=pattern->ptr[i],iptr<pattern->ptr[i+1],++i) {
          Degree[i] = Degree[i] + 1;
          diff = Label[i] - Label(pattern->index[iptr]);
          max_diff = MAX(max_diff,diff);
        }
        *profile + = max_diff;
        *bandwidth=MAX(*bandwidth,max_diff);
     }
}
/*  Finley_Reduce_findDiameter is the control procedure for finding the pseudo-diameter of pattern */
/*  as well as the level structure from each end                                                   */

/*  StartNode-        on input this is the node number of the first attempt at finding a diaMeter.  */
/*                    on output it  contains the actual number used.
/*  EndNode-          on output contains other end of diameter
/*  LevelToNode-       ARRAY CONTAINING LEVEL STRUCTURE WITH StartNode AS ROOT
/*  NodeToLevel-       ARRAY CONTAINING LEVEL STRUCTURE WITH EndNode AS ROOT
/*  IDFLT-       FLAG USED IN PICKING FINAL LEVEL STRUCTURE, SET
/*               =1 if WIDTH OF LevelToNode <= WIDTH OF NodeToLevel, OTHERWISE =2
/*  AssignedLevel,NodesInTree-     WORKING STORAGE
void Finley_Reduce_findDiameter(index_t *StartNode,index_t *EndNode,PPP pattern, NR, Degree, AssignedLevel, LevelToNode,NodeToLevel, NodesInTree, IDFLT)
      INTEGER NDSTK
      INTEGER FLAG, StartNode2, StartNode, EndNode
      COMMON /GRA/ N, IDPTH, ndeg_max
/* IT IS ASSUMED THAT THE LAST LEVEL HAS AT MOST 100 NODES.
      COMMON /CC/ NDLST(100)
      DIMENSION NDSTK(NR,1), Degree(1), AssignedLevel(1), LevelToNode(1), NodeToLevel(1),NodesInTree(1)
      FLAG = 0;
      MTW2 = N;
      StartNode2 = StartNode;
   10
     /* zero AssignedLevel to indicate all nodes are available to Finley_Reduce_dropTree */
     for (i=0,i<N,++i)
        AssignedLevel[i] = 0;

      available_NodesInTree = 1
     /* drop a tree from StartNode2 */
      CALL Finley_Reduce_dropTree(StartNode2, NDSTK, NR, AssignedLevel, NodesInTree, Degree, LastLevelWidth, BottomLevel,available_NodesInTree, max_LevelWidth, MTW2)
      if (FLAG<1) {
         FLAG = 1
   30    IDPTH = available_NodesInTree - 1
         MTW1 = max_LevelWidth
         /* copy level structure into LevelToNode */
         for (i=0,i<N,++i) 
            LevelToNode[i] = AssignedLevel[i];
         NDXN = 1;
         NDXL = 0;
         MTW2 = N;
         /* sort last level by degree  and store in ndlst
         CALL SORTDG(NDLST, NodesInTree(BottomLevel), NDXL, LastLevelWidth, Degree)
         StartNode2 = NDLST(1)
         GO TO 10
   } 
   50 if (IDPTH>=available_NodesInTree-1) GO TO 60
/* START AGAIN WITH NEW STARTING NODE
      StartNode = StartNode2
      GO TO 30
   60 if (max_LevelWidth>=MTW2) GO TO 80
      MTW2 = max_LevelWidth
      EndNode = StartNode2
/* STORE NARROWEST REVERSE LEVEL STRUCTURE IN NodeToLevel
      for (70 i=0,i<N,++i)
        NodeToLevel[i] = AssignedLevel[i]
   70 CONTINUE
   80 if (NDXN.EQ.NDXL) GO TO 90
/* TRY NEXT NODE IN NDLST
      NDXN = NDXN + 1
      StartNode2 = NDLST(NDXN)
      GO TO 10
   90 IDFLT = 1
      if (MTW2<=MTW1) IDFLT = 2
      RETURN
      END
*/
void Finley_Reduce_dropTree(index_t root, PPP pattern, NR, index_t *AssignedLevel,index_t *NodesInTree,
                            index_t *LastLevelWidth, index_t *BottomLevel,index_t *available_NodesInTree, index_t *max_LevelWidth, index_t max_LevelWidth_abort)

/*  Finley_Reduce_dropTree drops a tree in pattern from root */

/*  AssignedLevel- array of length length pattern->N indicating available nodes with zero entries. */
/*                 Finley_Reduce_dropTree enters level numbers assigned  during execution of this procedure */

/*  NodesInTree -  on output contains node numbers used in tree  (array of length pattern->N) */
/*                arranged by levels (NodesInTree[available_NodesInTree] contains root and NodesInTree(BottomLevel+LastLevelWidth-1) */
/*                contains last node entered)  */

/*  LastLevelWidth -    on output contains width of last level */
/*  BottomLevel    -    on output contains index into iwk of first node in last level
/*  max_LevelWidth -    on output contains the maximum level width
/*  available_NodesInTree-  on input the first available location in NodesInTree
/*                          usually one but if NodesInTree is used to store previous connected components, */
/*                          available_NodesInTree is next available location.  on output the total number of levels + 1 */
/*  max_LevelWidth_abort-       input param which triggers early return if max_LevelWidth becomes >= max_LevelWidth_abort */

  index_t ITOP = available_NodesInTree;
  index_t INOW = available_NodesInTree;

  *max_LevelWidth = 0;
  *BottomLevel = available_NodesInTree;
  *TopLevel = available_NodesInTree + 1;
  *available_NodesInTree = 1
  AssignedLevel(root) = 1
  NodesInTree(ITOP) = root
  while(max_LevelWidth<max_LevelWidth_abort && ITOP>=TopLevel) {
      available_NodesInTree++;
      while (INOW>=TopLevel) {
         NodesInTreeNOW = NodesInTree(INOW)
         for (j=pattern->iptr[NodesInTreeNOW],j<pattern->iptr[NodesInTreeNOW+1],++J) {
           ITEST = pattern->index[J];
           if (AssignedLevel(ITEST)==0) {
              AssignedLevel(ITEST) = available_NodesInTree;
              ITOP++;
              NodesInTree(ITOP) = ITEST
           }
         }
         INOW = INOW + 1
      }
      *LastLevelWidth = TopLevel - BottomLevel
      *max_LevelWidth = MAX(max_LevelWidth,LastLevelWidth);
      *BottomLevel = INOW
      *TopLevel = ITOP + 1
   }
}

      SUBROUTINE SORTDG(STK1, STK2, X1, X2, Degree)                       SOR   10
/* SORTDG SORTS STK2 BY DEGREE OF THE NODE AND ADDS IT TO THE END
/* OF STK1 IN ORDER OF LOWEST TO HIGHEST DEGREE.  X1 AND X2 ARE THE
/* NUMBER OF NODES IN STK1 AND STK2 RESPECTIVELY.
      INTEGER X1, X2, STK1, STK2, TEMP
      COMMON /GRA/ N, IDPTH, ndeg_max
      DIMENSION Degree(1), STK1(1), STK2(1)
      IND = X2
   10 ITEST = 0
      IND = IND - 1
      if (IND<1) GO TO 30
      for (20 I=1,IND
        J = I + 1
        ISTK2 = STK2[i]
        JSTK2 = STK2(J)
        if (Degree(ISTK2)<=Degree(JSTK2)) GO TO 20
        ITEST = 1
        TEMP = STK2[i]
        STK2[i] = STK2(J)
        STK2(J) = TEMP
   20 CONTINUE
      if (ITEST.EQ.1) GO TO 10
   30 for (40 I=1,X2
        X1 = X1 + 1
        STK1(X1) = STK2[i]
   40 CONTINUE
      RETURN
      END
      SUBROUTINE Finley_Reduce_setup(AssignedLevel, LevelToNode, NodeToLevel)                               SET   10
/* Finley_Reduce_setup COMPUTES THE REVERSE LEVELING INFO FROM NodeToLevel AND STORES
/* IT INTO NodeToLevel.  NACUM[i] IS INITIALIZED TO NODES/ITH LEVEL FOR NODES
/* ON THE PSEUDO-diameter OF THE GRAPH.  AssignedLevel IS INITIALIZED TO NON-
/* ZERO FOR NODES ON THE PSEUDO-diameter AND NODES IN A DifFERENT
/* COMPONENT OF THE GRAPH.
      COMMON /GRA/ N, IDPTH, ndeg_max
/* IT IS ASSUMED THAT THERE ARE AT MOST 100 LEVELS.
      COMMON /AssignedLevelW/ NHIGH(100), NLOW(100), NACUM(100)
      DIMENSION AssignedLevel(1), LevelToNode(1), NodeToLevel(1)
      for (10 I=1,IDPTH
        NACUM[i] = 0
   10 CONTINUE
      for (30 i=0,i<N,++i)
        AssignedLevel[i] = 1
        NodeToLevel[i] = IDPTH + 1 - NodeToLevel[i]
        ITEMP = NodeToLevel[i]
        if (ITEMP>IDPTH) GO TO 30
        if (ITEMP!=LevelToNode[i]) GO TO 20
        NACUM(ITEMP) = NACUM(ITEMP) + 1
        GO TO 30
   20   AssignedLevel[i] = 0
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION SORT2(DMY)                                       SOR   10
/* SORT2 SORTS SIZE AND STPT INTO DESCENDING ORDER ACCORDING TO
/* VALUES OF SIZE. ConnectedComponentCounter=NUMBER OF ENTRIES IN EACH ARRAY
      INTEGER TEMP, ConnectedComponents, SIZE, STPT, ConnectedComponentCounter
/* IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
      COMMON /CC/ ConnectedComponentCounter, SIZE(50), STPT(50)
      SORT2 = 0
      if (ConnectedComponentCounter.EQ.0) RETURN
      SORT2 = 1
      IND = ConnectedComponentCounter
   10 ITEST = 0
      IND = IND - 1
      if (IND<1) RETURN
      for (20 I=1,IND
        J = I + 1
        if (SIZE[i]>=SIZE(J)) GO TO 20
        ITEST = 1
        TEMP = SIZE[i]
        SIZE[i] = SIZE(J)
        SIZE(J) = TEMP
        TEMP = STPT[i]
        STPT[i] = STPT(J)
        STPT(J) = TEMP
   20 CONTINUE
      if (ITEST.EQ.1) GO TO 10
      RETURN
      END
      SUBROUTINE PIKAssignedLevel(LevelToNode, NodeToLevel, ConnectedComponents, IDFLT, DirectionLargestComponent)             PIK   10
/* PIKAssignedLevel CHOOSES THE LEVEL STRUCTURE  USED IN NUMBERING GRAPH
/* LevelToNode-       ON INPUT CONTAINS FORWARD LEVELING INFO
/* NodeToLevel-       ON INPUT CONTAINS REVERSE LEVELING INFO
/*              ON OUTPUT THE FINAL LEVEL STRUCTURE CHOSEN
/* ConnectedComponents-      ON INPUT CONTAINS CONNECTED COMPONENT INFO
/* IDFLT-       ON INPUT =1 if WDTH LevelToNode<=WDTH NodeToLevel, =2 OTHERWISE
/* NHIGH        KEEPS TRACK OF LEVEL WIDTHS FOR HIGH NUMBERING
/* NLOW-        KEEPS TRACK OF LEVEL WIDTHS FOR LOW NUMBERING
/* NACUM-       KEEPS TRACK OF LEVEL WIDTHS FOR CHOSEN LEVEL STRUCTURE
/* ConnectedComponentCounter-          NUMBER OF CONNECTED COMPONENTS
/* SIZE[i]-     SIZE OF ITH CONNECTED COMPONENT
/* STPT[i]-     INDEX INTO ConnectedComponentsE OF 1ST NODE IN ITH CON COMPT
/* DirectionLargestComponent-       FLAG WHICH INDICATES WHICH WAY THE LARGEST CONNECTED
/*              COMPONENT FELL.  =+1 if LOW AND -1 if HIGH
      INTEGER ConnectedComponents, SIZE, STPT, ConnectedComponentCounter, END
      COMMON /GRA/ N, IDPTH, ndeg_max
/* IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 COMPONENTS AND
/* THAT THERE ARE AT MOST 100 LEVELS.
      COMMON /AssignedLevelW/ NHIGH(100), NLOW(100), NACUM(100)
      COMMON /CC/ ConnectedComponentCounter, SIZE(50), STPT(50)
      DIMENSION LevelToNode(1), NodeToLevel(1), ConnectedComponents(1)
/* FOR EACH CONNECTED COMPONENT DO
      for (80 I=1,ConnectedComponentCounter
        J = STPT[i]
        END = SIZE[i] + J - 1
/* SET NHIGH AND NLOW EQUAL TO NACUM
        for (10 K=1,IDPTH
          NHIGH(K) = NACUM(K)
          NLOW(K) = NACUM(K)
   10   CONTINUE
/* UPDATE NHIGH AND NLOW FOR EACH NODE IN CONNECTED COMPONENT
        for (20 K=J,END
          INODE = ConnectedComponents(K)
          available_NodesInTreeH = LevelToNode(INODE)
          NHIGH(available_NodesInTreeH) = NHIGH(available_NodesInTreeH) + 1
          available_NodesInTreeL = NodeToLevel(INODE)
          NLOW(available_NodesInTreeL) = NLOW(available_NodesInTreeL) + 1
   20   CONTINUE
        MAX1 = 0
        MAX2 = 0
/* SET MAX1=LARGEST NEW NUMBER IN NHIGH
/* SET MAX2=LARGEST NEW NUMBER IN NLOW
        for (30 K=1,IDPTH
          if (2*NACUM(K).EQ.NLOW(K)+NHIGH(K)) GO TO 30
          if (NHIGH(K)>MAX1) MAX1 = NHIGH(K)
          if (NLOW(K)>MAX2) MAX2 = NLOW(K)
   30   CONTINUE
/* SET IT= NUMBER OF LEVEL STRUCTURE TO BE USED
        IT = 1
        if (MAX1>MAX2) IT = 2
        if (MAX1.EQ.MAX2) IT = IDFLT
        if (IT.EQ.2) GO TO 60
        if (I.EQ.1) DirectionLargestComponent = -1
/* COPY LevelToNode INTO NodeToLevel FOR EACH NODE IN CONNECTED COMPONENT
        for (40 K=J,END
          INODE = ConnectedComponents(K)
          NodeToLevel(INODE) = LevelToNode(INODE)
   40   CONTINUE
/* UPDATE NACUM TO BE THE SAME AS NHIGH
        for (50 K=1,IDPTH
          NACUM(K) = NHIGH(K)
   50   CONTINUE
        GO TO 80
/* UPDATE NACUM TO BE THE SAME AS NLOW
   60   for (70 K=1,IDPTH
          NACUM(K) = NLOW(K)
   70   CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE NUMBER(StartNode2, NUM, NDSTK, NodeToLevel, Degree, NewLabel, AssignedLevelST,     NUM   10
     * LSTPT, NR, Flag, newBandwidth, newProfile, IPFA, DirectionLargestComponent)
/*  NUMBER PRODUCES THE NUMBERING OF THE GRAPH FOR MIN BANDWIDTH
/*  StartNode2-         ON INPUT THE NODE TO BEGIN NUMBERING ON
/*  NUM-         ON INPUT AND OUTPUT, THE NEXT AVAILABLE NUMBER
/*  NodeToLevel-       THE LEVEL STRUCTURE TO BE USED IN NUMBERING
/*  NewLabel-       THE ARRAY USED TO STORE THE NEW NUMBERING
/*  AssignedLevelST-       ON OUTPUT CONTAINS LEVEL STRUCTURE
/*  LSTPT[i]-    ON OUTPUT, INDEX INTO AssignedLevelST TO FIRST NODE IN ITH AssignedLevel
/*               LSTPT(I+1) - LSTPT[i] = NUMBER OF NODES IN ITH AssignedLevel
/*  Flag-        =+1 if StartNode2 IS FORWARD END OF PSEUDO-diameter
/*               =-1 if StartNode2 IS REVERSE END OF PSEUDO-diameter
/*  newBandwidth-        BANDWIDTH OF NEW NUMBERING COMPUTED BY NUMBER
/*  newProfile-        PROFILE OF NEW NUMBERING COMPUTED BY NUMBER
/*  IPFA-        WORKING STORAGE USED TO COMPUTE PROFILE AND BANDWIDTH
/*  DirectionLargestComponent-       INDICATES STEP DIRECTION USED IN NUMBERING(+1 OR -1)
/* USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
      INTEGER NDSTK
      INTEGER StartNode2, STKA, STKB, STKC, STKD, XA, XB, ConnectedComponentCounter, XD, CX, END,
     * NewLabel, TEST
      COMMON /GRA/ N, IDPTH, ndeg_max
/* THE STORAGE IN COMMON BLOCKS CC AND AssignedLevelW IS NOW FREE AND CAN
/* BE USED FOR STACKS.
      COMMON /AssignedLevelW/ STKA(100), STKB(100), STKC(100)
      COMMON /CC/ STKD(100)
      DIMENSION IPFA(1)
      DIMENSION NDSTK(NR,1), NodeToLevel(1), Degree(1), NewLabel(1), AssignedLevelST(1),
     * LSTPT(1)
/* SET UP AssignedLevelST AND LSTPT FROM NodeToLevel
      for (10 i=0,i<N,++i)
        IPFA[i] = 0
   10 CONTINUE
      NSTPT = 1
      for (30 I=1,IDPTH
        LSTPT[i] = NSTPT
        for (20 J=1,N
          if (NodeToLevel(J)!=I) GO TO 20
          AssignedLevelST(NSTPT) = J
          NSTPT = NSTPT + 1
   20   CONTINUE
   30 CONTINUE
      LSTPT(IDPTH+1) = NSTPT
/* STKA, STKB, STKC AND STKD ARE STACKS WITH POINTERS
/* XA,XB,ConnectedComponentCounter, AND XD.  CX IS A SPECIAL POINTER INTO STKC WHICH
/* INDICATES THE PARTICULAR NODE BEING PROCESSED.
/* available_NodesInTree KEEPS TRACK OF THE LEVEL WE ARE WORKING AT.
/* INITIALLY STKC CONTAINS ONLY THE INITIAL NODE, StartNode2.
      available_NodesInTree = 0
      if (Flag<0) available_NodesInTree = IDPTH + 1
      ConnectedComponentCounter = 1
      STKC(ConnectedComponentCounter) = StartNode2
   40 CX = 1
      XD = 0
      available_NodesInTree = available_NodesInTree + Flag
      LST = LSTPT(available_NodesInTree)
      LND = LSTPT(available_NodesInTree+1) - 1
/* BEGIN PROCESSING NODE STKC(CX)
   50 IPRO = STKC(CX)
      NewLabel(IPRO) = NUM
      NUM = NUM + DirectionLargestComponent
      END = Degree(IPRO)
      XA = 0
      XB = 0
/* CHECK ALL ADJACENT NODES
      for (80 I=1,END
        TEST = NDSTK(IPRO,I)
        INX = NewLabel(TEST)
/* ONLY NODES NOT NUMBERED OR ALREADY ON A STACK ARE ADDED
        if (INX.EQ.0) GO TO 60
        if (INX<0) GO TO 80
/* for (PRELIMINARY BANDWIDTH AND PROFILE CALCULATIONS
        NBW = (NewLabel(IPRO)-INX)*DirectionLargestComponent
        if (DirectionLargestComponent>0) INX = NewLabel(IPRO)
        if (IPFA(INX)<NBW) IPFA(INX) = NBW
        GO TO 80
   60   NewLabel(TEST) = -1
/* PUT NODES ON SAME LEVEL ON STKA, ALL OTHERS ON STKB
        if (NodeToLevel(TEST).EQ.NodeToLevel(IPRO)) GO TO 70
        XB = XB + 1
        STKB(XB) = TEST
        GO TO 80
   70   XA = XA + 1
        STKA(XA) = TEST
   80 CONTINUE
/* SORT STKA AND STKB INTO INCREASING DEGREE AND ADD STKA TO STKC
/* AND STKB TO STKD
      if (XA.EQ.0) GO TO 100
      if (XA.EQ.1) GO TO 90
      CALL SORTDG(STKC, STKA, ConnectedComponentCounter, XA, Degree)
      GO TO 100
   90 ConnectedComponentCounter = ConnectedComponentCounter + 1
      STKC(ConnectedComponentCounter) = STKA(XA)
  100 if (XB.EQ.0) GO TO 120
      if (XB.EQ.1) GO TO 110
      CALL SORTDG(STKD, STKB, XD, XB, Degree)
      GO TO 120
  110 XD = XD + 1
      STKD(XD) = STKB(XB)
/* BE SURE TO PROCESS ALL NODES IN STKC
  120 CX = CX + 1
      if (ConnectedComponentCounter>=CX) GO TO 50
/* WHEN STKC IS EXHAUSTED LOOK FOR MIN DEGREE NODE IN SAME LEVEL
/* WHICH HAS NOT BEEN PROCESSED
      MAX = ndeg_max + 1
      StartNode2 = N + 1
      for (130 I=LST,LND
        TEST = AssignedLevelST[i]
        if (NewLabel(TEST)!=0) GO TO 130
        if (Degree(TEST)>=MAX) GO TO 130
        NewLabel(StartNode2) = 0
        NewLabel(TEST) = -1
        MAX = Degree(TEST)
        StartNode2 = TEST
  130 CONTINUE
      if (StartNode2.EQ.N+1) GO TO 140
      ConnectedComponentCounter = ConnectedComponentCounter + 1
      STKC(ConnectedComponentCounter) = StartNode2
      GO TO 50
/* if STKD IS EMPTY WE ARE DONE, OTHERWISE COPY STKD ONTO STKC
/* AND BEGIN PROCESSING NEW STKC
  140 if (XD.EQ.0) GO TO 160
      for (150 I=1,XD
        STKC[i] = STKD[i]
  150 CONTINUE
      ConnectedComponentCounter = XD
      GO TO 40
/* for (FINAL BANDWIDTH AND PROFILE CALCULATIONS
  160 for (170 i=0,i<N,++i)
        if (IPFA[i]>newBandwidth) newBandwidth = IPFA[i]
        newProfile = newProfile + IPFA[i]
  170 CONTINUE
      RETURN
      END
