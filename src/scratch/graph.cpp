#include "graph.h"
#include <ctime>
using namespace Eigen;

Graph::Graph(int V){
    this->V = V;
    adj = new list<int>[V];
    Nedges=0;}

void Graph::addEdge(int v, int w){
    adj[v].push_back(w);    // adjacency list
    edge.push_back(graphedge(v,w));    // edge vector
    Nedges++;
    Nnodes=MAX(Nnodes,w++);
    Nnodes=MAX(Nnodes,v++);}   // constraints structure

void Graph::addConstraint(int v, int u, vector<double> QQ, vector<double> t){
    vector<int> I; I.push_back(v); I.push_back(u);
    vector<double> e; e.assign(4,0);
    //e.push_back(0);e.push_back(0);e.push_back(0);e.push_back(0);
    C.push_back(Constraints(I,QQ,e,t));
    Nedges++;
    Nnodes=MAX(Nnodes,u++);
    Nnodes=MAX(Nnodes,v++);}   // constraints structure

void Graph::initlialisePose( ){
    vector<double> Q; Q.assign(4,0); Q[0] = 1;
    vector<double> t; t.assign(3,0);
    //Q.push_back(1);Q.push_back(0);Q.push_back(0);Q.push_back(0);
    //t.push_back(0);t.push_back(0);t.push_back(0);
    // initialise M
    //std::cout << nodes << std::endl;
    for(int j=0;j<Nnodes;j++){ M.push_back(Pose(Q,t)); }}
/* to update Q or t:
 * for(int j=0;j<V;j++){
 *      Pose &pose = M[j]; // grab a reference to one of the poses
 *      pose.Q[0] = 0; // updates the first valie of Q
 *      pose.t[0] = 0;} // updates the first valie of t
 */

// The function to do DFS traversal. It uses SCCUtil()
// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
//             discovery time) that can be reached from subtree
//             rooted with current vertex
// *st -- >> To store all the connected ancestors (could be part
//           of SCC)
// stackMember[] --> bit/index array for faster check whether
//                  a node is in stack
void Graph::SCC(){
    int *disc = new int[V];
    int *low = new int[V];
    bool *stackMember = new bool[V];
    stack<int> *st = new stack<int>();
    // Initialize disc and low, and stackMember arrays
    for (int i = 0; i<V; i++){ disc[i] = NIL; low[i] = NIL;
    stackMember[i] = false;}
    // Call the recursive helper function to find strongly
    // connected components in DFS tree with vertex 'i'
    for (int i=0;i<V;i++){
        if (disc[i]==NIL){SCCUtil(i,disc,low,st,stackMember);}}
    // fix indecies to count incrementally get outputs similar to matlab's
    // graphconncomp.m (S and CC)
    int temp=sub[0].index,index=0;
    int m=0,j=0,k=0;
    int i=0; CC.assign(Nnodes,-1);
    while (i<sub.size()){
        if (temp==sub[i].index){k++;
        if (k>m){m=k;SS=sub[i].index;//or (index)
        }}
        else {k=1;temp=sub[i].index;index++;}
        CC[sub[i].vertex]=sub[i].index;//or (index)
        i++;}
    // finding the edges with nodes in the SCC
    //for (i=0;i<Nedges;i++){
    //    list<int>::iterator it;
    //    for (it = adj[i].begin(); it != adj[i].end(); ++it){
    //        if (CC[*it]==SS&CC[i]==SS){JJ[i]=1;k++;}
    //        else {JJ[k]=0;k++;}}}
    for (i=0;i<Nedges;i++){
        if (CC[edge[i].i]==SS&CC[edge[i].j]==SS){JJ.push_back(1);}
        else {JJ.push_back(0);}}}

void Graph::SCCUtil(int u, int disc[], int low[], stack<int> *st,
        bool stackMember[]){
    // A static variable is used for simplicity, we can avoid use
    // of static variable by passing a pointer.
    static int time = 0;
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    st->push(u);
    stackMember[u] = true;
    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i){
        int v = *i;  // v is current adjacent of 'u'
        // If v is not visited yet, then recur for it
        if (disc[v] == -1){
            SCCUtil(v, disc, low, st, stackMember);
            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 (per above discussion on Disc and Low value)
            low[u]  = min(low[u], low[v]);}
        // Update low value of 'u' only of 'v' is still in stack
        // (i.e. it's a back edge, not cross edge).
        // Case 2 (per above discussion on Disc and Low value)
        else if (stackMember[v]==true){low[u]  = min(low[u], disc[v]);}}
    // head node found, pop the stack and print an SCC
    int w = 0;  // To store stack extracted vertices
    if (low[u] == disc[u]){
        while (st->top() != u){
            w = (int) st->top();
            //cout << w << " ";
            sub.push_back(subgraph(w,u));
            stackMember[w] = false;
            st->pop();}
        w = (int) st->top();
        //cout << w << "\n";
        sub.push_back(subgraph(w,u));
        stackMember[w] = false;
        st->pop();}}

// vector<double> Graph::RtoQuaternion(vector<double> R){
//     Q(1)=(R.coeffRef(1,1)+R.coeffRef(2,2)+R.coeffRef(3,3)-1)/2;
//     Q(2)=(R.coeffRef(3,2)-R.coeffRef(2,3))/2;
//     Q(3)=(R.coeffRef(1,3)-R.coeffRef(3,1))/2;
//     Q(4)=(R.coeffRef(2,1)-R.coeffRef(1,2))/2;
//     Q(1)=sqrt((Q(1)+1)/2);
//     Q(2)=(Q(2)/Q(1))/2;
//     Q(3)=(Q(3)/Q(1))/2;
//     Q(4)=(Q(4)/Q(1))/2;}

void Graph::QuaternionInit( ){
    // initialise M.Q from a spanning-tree
    vector<int> it;
    it.assign(Nnodes,0); // fill the vector i with all zeros
    int k,i,j,a=0,SpanFlag,num_nodes=1;
    it[a] = 1; // first node
    //num_nodes = std::accumulate(it.begin(),it.end(),0);
    while(num_nodes<Nnodes){
        SpanFlag=0;
        for(k=0;k<C.size();k++){
            i = C[k].I[0]; j = C[k].I[1];
            if (it[i]==1&&it[j]==0){it[j] = 1; SpanFlag=1;
            M[j].Q=QuaternionFwd(M[i].Q,C[k].QQ);}
            if (it[i]==0&&it[j]==1){it[i] = 1; SpanFlag=1;
            M[i].Q=QuaternionRvs(M[j].Q,C[k].QQ);}}
        num_nodes = std::accumulate(it.begin(),it.end(),0);
        if (SpanFlag==0&&num_nodes<Nnodes){
            std::cout << "Relative rotations DO NOT SPAN all the nodes in the VIEW GRAPH" << std::endl;
            std::cout << "Number of nodes in Spanning Tree = " << num_nodes << std::endl;
            std::cout << "Connected Nodes are : ";
            for (k=0;k<it.size();k++){ std::cout << it[k] << " ";} std::cout << std::endl;
            std::cout << "Remove extra nodes and retry" << std::endl;
            exit(1);}}}

void Graph::QuaternionBoxMedian(){
    
    std::clock_t start1 = std::clock();
    double duration;
    
    /* Initialise quaternions */
    //QuaternionInit( );
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"[RotAvg] initialisation : "<< duration <<"s.\n";
    
    /* Form system matrix (A) and (AtA) */
    int i,j;
    int n=Nnodes-1,m=Nedges;
    vector<T> AList,AtAList;
    SpMat A(m,n),AtA(n*n,m);
    for(int k=0;k<C.size();k++){
        i = C[k].I[0]; j = C[k].I[1];
        /* The first node orientation is fixed at Q1=[1;0;0;0], hence the
         * first node is also removed from the system. The indecies become
         * (i -> i-1) and (j -> j-1)*/
        if (i>0){AList.push_back(T(k,i-1,-1));
        AtAList.push_back(T(n*(i-1)+i-1,k, 1));
        if (j>0){AtAList.push_back(T(n*(i-1)+j-1,k,-1));}}
        if (j>0){AList.push_back(T(k,j-1, 1));
        AtAList.push_back(T(n*(j-1)+j-1,k, 1));
        if (i>0){AtAList.push_back(T(n*(j-1)+i-1,k,-1));}}}
    A.setFromTriplets(AList.begin(),AList.end());
    AtA.setFromTriplets(AtAList.begin(),AtAList.end());
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"[RotAvg] matrix AtA : \t\t"<< duration <<"s.\n";
    
    /* Notice that during minimisation the system matrix A is not changing.
     * A better implementation will consider removing and adding
     * constraints to increase the robustness of the solution against
     * outlier motions */
    int Iteration=0, L1Step=2, maxIters=100;
    double changeThreshold=0.001, temp, t;
    double s2, theta, score=1000;
    vector<double> e; e.assign(4,0);
    vector<double> QQ; QQ.assign(4,0);
    
    //std::cout << "Iteration" << "\t" << "score" << std::endl;
    
   while(((score>=changeThreshold)||(L1Step<2))&&(Iteration<maxIters)){
        
        if(score<changeThreshold){
            L1Step=L1Step*4;
            changeThreshold=changeThreshold/100;}
        
        /* Form rotation error matrix (B) and solve the system */
        vector<T> BList; SpMat B(m,3);
        for(int k=0;k<C.size();k++){
            i = C[k].I[0]; j = C[k].I[1];
            e=QuaternionFwd(M[i].Q, C[k].QQ); // e=Qij*Qi
            e=QuaternionRvs(e     , M[j].Q ); // e=inv(Qj)*e=inv(Qj)*Qij*Qi
            s2=sqrt(e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);
            e[0]=2*atan2(s2,e[0]);
            if (e[0]<-pi){e[0]=e[0]+2*pi;}
            if (e[0]>=pi){e[0]=e[0]-2*pi;}
            t = e[1]*e[0]/s2; if (!isnan(t)){BList.push_back(T(k,0,t));}
            t = e[2]*e[0]/s2; if (!isnan(t)){BList.push_back(T(k,1,t));}
            t = e[3]*e[0]/s2; if (!isnan(t)){BList.push_back(T(k,2,t));}}
        B.setFromTriplets(BList.begin(), BList.end());
        
        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        //std::cout<<"[RotAvg] matrix B : \t\t"<< duration <<"s.\n";
        
        /* Decoding via linear programming: Solve min_x(||b-Ax||_1) */
        MatXd W  = MatXd::Zero(n,3);
        MatXd Bd = MatXd(B);
        W.col(0) = l1decode_pd(W.col(0), A, Bd.col(0), AtA, L1Step);
        W.col(1) = l1decode_pd(W.col(1), A, Bd.col(1), AtA, L1Step);
        W.col(2) = l1decode_pd(W.col(2), A, Bd.col(2), AtA, L1Step);
        
        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        //std::cout<<"[RotAvg] Median solving : \t"<< duration <<"s.\n";
        
        /* Augment the first node */
        MatXd W2 = MatXd::Zero(n+1,4);
        W2.bottomRightCorner(n, 3) = W;
        W2(0,0) = 1;
        /* score */
        score=0;
        for (int i=0;i<W2.innerSize();i++){
            temp=0;
            for (int j=1;j<W2.outerSize();j++){
                temp+=W2.coeffRef(i,j)*W2.coeffRef(i,j);}
            temp=sqrt(temp);
            if (i>0){ score=MAX(score,temp); }// score without the first camera
            for (int j=0;j<W2.outerSize();j++){
                if (j==0){ W2.coeffRef(i,j)=cos(temp/2); }
                else { W2.coeffRef(i,j)=W2.coeffRef(i,j)*sin(temp/2)/temp; }
                if (isnan(W2.coeffRef(i,j))){ W2.coeffRef(i,j)=0;}
                QQ[j]=W2.coeffRef(i,j);} // error Quaternion
            M[i].Q=QuaternionFwd(QQ,M[i].Q);} // update Quaternion
        
        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        //std::cout<<"[RotAvg] Update : \t\t"<< duration <<"s.\n";
        
        //std::cout << ++Iteration << "\t" << score << std::endl;}}
       ++Iteration;}
    
    duration = ( std::clock() - start1 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"[RotAvg] Median time : \t\t"<< duration <<"s.\n";
    
}

void Graph::QuaternionRobustMean(){
    
    std::clock_t start2 = std::clock();
    double duration;
    
    /* Form system matrix (A) and (AtA) */
    int i,j;
    int n=Nnodes,m=Nedges;
    
    /* Creat A as dense matrix */
//     MatXd Ad = MatXd::Zero(m,n-1);
//     for(int k=0;k<C.size();k++){
//         i = C[k].I[0]; j = C[k].I[1];
//         /* The first node orientation is fixed at Q1=[1;0;0;0], hence the
//          * first node is also removed from the system. The indecies become
//          * (i -> i-1) and (j -> j-1)*/
//         if (i>0){Ad.coeffRef(k,i-1)=-1;}
//         if (j>0){Ad.coeffRef(k,j-1)= 1;}}
    
    /* Creat A as sparse matrix */
    SpMat As(m,n-1);
    vector<T> AList;
    for(int k=0;k<C.size();k++){
        i = C[k].I[0]; j = C[k].I[1];
        if (i>0){AList.push_back(T(k,i-1,-1));}
        if (j>0){AList.push_back(T(k,j-1, 1));}}
    As.setFromTriplets(AList.begin(),AList.end());
    //MatXd Ad = MatrixXd(As);
    
    int Iteration=0, maxIters=100;
    double SIGMA=5*pi/180; // Degrees to radians
    double changeThreshold=0.001, temp, t;
    double s2, theta, score=1000;
    vector<double> e;   e.assign(4,0);
    vector<double> QQ; QQ.assign(4,0);
    VecXd Weights = MatXd::Ones(m,1);
    
    //std::cout << "Iteration" << "\t" << "score" << std::endl;
    
    while((score>=changeThreshold)&&(Iteration<maxIters)){
        
        /* Form rotation error matrix (B) and solve the system */
        
        /* as Dense */
        MatXd B = MatXd::Zero(m,3);
        MatXd B2 = MatXd::Zero(m,3);
        for(int k=0;k<C.size();k++){
            i = C[k].I[0]; j = C[k].I[1];
            e=QuaternionFwd(M[i].Q, C[k].QQ); // e=Qij*Qi
            e=QuaternionRvs(e     , M[j].Q ); // e=inv(Qj)*e=inv(Qj)*Qij*Qi
            s2=sqrt(e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);
            e[0]=2*atan2(s2,e[0]);
            if (e[0]<-pi){e[0]=e[0]+2*pi;}
            if (e[0]>=pi){e[0]=e[0]-2*pi;}
            t = e[1]*e[0]/s2;
            if (!isnan(t)){B.coeffRef(k,0)=t;
            B2.coeffRef(k,0)=Weights.coeffRef(k)*t;}
            t = e[2]*e[0]/s2;
            if (!isnan(t)){B.coeffRef(k,1)=t;
            B2.coeffRef(k,1)=Weights.coeffRef(k)*t;}
            t = e[3]*e[0]/s2;
            if (!isnan(t)){B.coeffRef(k,2)=t;
            B2.coeffRef(k,2)=Weights.coeffRef(k)*t;}}

        /* as Sparse */
//         vector<T> BList, B2List;
//         SpMat Bs(m,3), B2s(m,3);
//         for(int k=0;k<C.size();k++){
//             i = C[k].I[0]; j = C[k].I[1];
//             e=QuaternionFwd(M[i].Q, C[k].QQ); // e=Qij*Qi
//             e=QuaternionRvs(e     , M[j].Q ); // e=inv(Qj)*e=inv(Qj)*Qij*Qi
//             s2=sqrt(e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);
//             e[0]=2*atan2(s2,e[0]);
//             if (e[0]<-pi){e[0]=e[0]+2*pi;}
//             if (e[0]>=pi){e[0]=e[0]-2*pi;}
//             t = e[1]*e[0]/s2; if (!isnan(t)){
//                 BList.push_back(T(k,0,t));
//                 B2List.push_back(T(k,0,Weights.coeffRef(k)*t));}
//             t = e[2]*e[0]/s2; if (!isnan(t)){
//                 BList.push_back(T(k,1,t));
//                 B2List.push_back(T(k,1,Weights.coeffRef(k)*t));}
//             t = e[3]*e[0]/s2; if (!isnan(t)){
//                 BList.push_back(T(k,2,t));
//                 B2List.push_back(T(k,2,Weights.coeffRef(k)*t));}}
//         Bs.setFromTriplets(BList.begin(), BList.end());
//         B2s.setFromTriplets(B2List.begin(), B2List.end());
//         MatXd B  = MatrixXd(Bs);
//         MatXd B2 = MatrixXd(B2s);
        
        /* Weights diagonal matrix */
        /* using MatXd --------------------------------------- speed -- */
//         MatXd Ws = MatXd::Zero(m,m);
//         for (int k=0; k<m; k++){
//             Ws.coeffRef(k,k)=Weights.coeffRef(k);}
        /* using asDiagonal ---------------------------------- speed -  */
//         MatXd Wd = Weights.asDiagonal();
//         SpMat Ws = Wd.sparseView();
        /* filling a sparse matrix --------------------------- speed ++ */
        SpMat Ws(m,m);
        vector<T> WList;
        for(int k=0;k<m;k++){
            WList.push_back(T(k,k,Weights.coeffRef(k)));}
        Ws.setFromTriplets(WList.begin(),WList.end());
        
        /* solve */
        /* multiplication using MatXd ------------------------ speed - */
//        MatXd A2d = Wd*Ad;
        /* multiplication using SpMat ------------------------ speed + */
        SpMat A2s = Ws*As;
        MatXd A2d = MatrixXd(A2s);
        // x = A \ b.
        //MatXd W = A2.fullPivLu().solve(B2);          // speed --
        //MatXd W = A2.colPivHouseholderQr().solve(B2);// speed -
        MatXd W = A2d.householderQr().solve(B2);       // speed ++
        
        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        //std::cout<<"[RotAvg] Robust solving : \t"<< duration <<"s.\n";
        
        // W2
        MatXd W2 = MatXd::Zero(n, 4);
        W2.bottomRightCorner(n-1, 3) = W;
        W2(0,0) = 1;
        
        // more like residuals
        MatXd Ad = MatrixXd(As);
        MatXd E = Ad*W - B;
        for (int i=0; i<E.innerSize(); i++){
            temp=0;
            for (int j=0; j<E.outerSize(); j++){
                temp+=E.coeffRef(i,j)*E.coeffRef(i,j);}
            Weights.coeffRef(i)=SIGMA/(temp+SIGMA*SIGMA);}
        
        /* score */
        score=0;
        for (int i=0;i<W2.innerSize();i++){
            temp=0;
            for (int j=1;j<W2.outerSize();j++){
                temp+=W2.coeffRef(i,j)*W2.coeffRef(i,j);}
            temp=sqrt(temp);
            if (i>0){ score+=temp/n;}// score without the first camera
            for (int j=0;j<W2.outerSize();j++){
                if (j==0){ W2.coeffRef(i,j)=cos(temp/2); }
                else { W2.coeffRef(i,j)=W2.coeffRef(i,j)*sin(temp/2)/temp; }
                if (isnan(W2.coeffRef(i,j))){ W2.coeffRef(i,j)=0;}
                QQ[j]=W2.coeffRef(i,j);} // error Quaternion
            M[i].Q=QuaternionFwd(QQ,M[i].Q);} // update Quaternion
        //std::cout << ++Iteration << "\t" << score << std::endl;}}
        ++Iteration;}
    
    duration = ( std::clock() - start2 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"[RotAvg] Robust time : \t\t"<< duration <<"s.\n";}

vector<double> Graph::QuaternionFwd(vector<double> Qi, vector<double> Qij){
    // Qj=Qij*Qi
    vector<double> Qj;
    Qj.push_back(Qij[0]*Qi[0]-Qij[1]*Qi[1]-Qij[2]*Qi[2]-Qij[3]*Qi[3]);
    Qj.push_back(Qij[1]*Qi[0]+Qij[0]*Qi[1]-Qij[3]*Qi[2]+Qij[2]*Qi[3]);
    Qj.push_back(Qij[2]*Qi[0]+Qij[3]*Qi[1]+Qij[0]*Qi[2]-Qij[1]*Qi[3]);
    Qj.push_back(Qij[3]*Qi[0]-Qij[2]*Qi[1]+Qij[1]*Qi[2]+Qij[0]*Qi[3]);
    return Qj;}

vector<double> Graph::QuaternionRvs(vector<double> Qj, vector<double> Qij){
    // Qi=Qij*Qj
    vector<double> Qi;
    Qi.push_back(-Qij[0]*Qj[0]-Qij[1]*Qj[1]-Qij[2]*Qj[2]-Qij[3]*Qj[3]);
    Qi.push_back( Qij[1]*Qj[0]-Qij[0]*Qj[1]-Qij[3]*Qj[2]+Qij[2]*Qj[3]);
    Qi.push_back( Qij[2]*Qj[0]+Qij[3]*Qj[1]-Qij[0]*Qj[2]-Qij[1]*Qj[3]);
    Qi.push_back( Qij[3]*Qj[0]-Qij[2]*Qj[1]+Qij[1]*Qj[2]-Qij[0]*Qj[3]);
    return Qi;}

/*************************************************************************
 * l1decode_pd.m
 *
 * Decoding via linear programming.
 * Solve
 * min_x  ||b-Ax||_1 .
 *
 * Recast as the linear program
 * min_{x,u} sum(u)  s.t.  -Ax - u + y <= 0
 *                          Ax - u - y <= 0
 * and solve using primal-dual interior point method.
 *
 * Usage: xp = l1decode_pd(x0, A, At, y, pdtol, pdmaxiter, cgtol, cgmaxiter)
 *
 * x0 - Nx1 vector, initial point.
 *
 * A -  Either a handle to a function that takes a N vector and returns a M
 *      vector, or a MxN matrix.  If A is a function handle, the algorithm
 *      operates in "largescale" mode, solving the Newton systems via the
 *      Conjugate Gradients algorithm.
 *
 * At - Handle to a function that takes an M vector and returns an N vector.
 *      If A is a matrix, At is ignored.
 *
 * y  - Mx1 observed code (M > N).
 *
 * pdtol -  Tolerance for primal-dual algorithm (algorithm terminates if
 *          the duality gap is less than pdtol).
 *          Default = 1e-3.
 *
 * pdmaxiter -  Maximum number of primal-dual iterations.
 *              Default = 50
 *************************************************************************/
MatXd Graph::l1decode_pd(VecXd x0, SpMat A, VecXd y, SpMat AtA, int pdmaxiter){
    double pdtol  = 1e-3;
    double alpha  = 0.01;
    double beta   = 0.5;
    double mu     = 10;
    int n = Nnodes-1, m = Nedges;
    
    /* u */
    VecXd x  = x0;     // (nodes-1) x 1
    //MatXd Ad = MatrixXd(A);
    VecXd Ax = A*x0;   // V x 1
    VecXd u  = y-Ax;
    u = u.cwiseAbs();
    u = 0.95*u + 0.1*u.maxCoeff()*MatXd::Ones(m,1);
    
    /* tau */
    vector<double> fu1, fu2, lamu1, lamu2; // V x 1
    double sdg=0;
    for (int k=0; k<m; k++){
        fu1.push_back( Ax.coeffRef(k) - y.coeffRef(k) - u.coeffRef(k));
        fu2.push_back(-Ax.coeffRef(k) + y.coeffRef(k) - u.coeffRef(k));
        lamu1.push_back(-1/fu1[k]);
        lamu2.push_back(-1/fu2[k]);
        sdg = sdg - fu1[k]*lamu1[k] - fu2[k]*lamu2[k];}
    double tau = mu*2*m/sdg;
    
    /* resnorm */
    vector<double> Atv; // (nodes-1) x 1
    double resnorm=0;
    for (int k=0; k<n; k++){
        Atv.push_back(0);
        for (SparseMatrix<double,ColMajor>::InnerIterator it(A,k);it;++it){
            Atv[k] += A.coeffRef(it.row(),it.col())*(lamu1[it.row()]-lamu2[it.row()]);}
        resnorm += Atv[k]*Atv[k];}
    double temp;
    for (int k=0; k<m; k++){
        temp    = -lamu1[k]*fu1[k]-(1/tau);
        resnorm += temp*temp;
        temp    = -lamu2[k]*fu2[k]-(1/tau);
        resnorm += temp*temp;
        temp    = -lamu1[k]-lamu2[k]+1;
        resnorm += temp*temp;}
    resnorm = sqrt(resnorm);
    
    /* While loop minimising resnorm */
    int pditer = 0;
    bool done = (sdg < pdtol) | (pditer >= pdmaxiter);
    
    vector<double> w1,w2;
    vector<double> sig1,sig2,sigx,sigy;
    vector<double> Atdv,du,dlamu1,dlamu2;
    vector<double> Atvp,fu1p,fu2p,lamu1p,lamu2p;
    VecXd w1p(n),xp(n),up(m),Axp(m);
    //SpMat H11p(n,n);
    
    while (!done){
        pditer = pditer + 1;
        
        /* w1p */
        for (int k=0; k<m; k++){
            w2.push_back(-1-(1/tau)*(1/fu1[k]+1/fu2[k]));
            sig1.push_back(-lamu1[k]/fu1[k]-lamu2[k]/fu2[k]);
            sig2.push_back(lamu1[k]/fu1[k]-lamu2[k]/fu2[k]);
            sigx.push_back(sig1[k]-(sig2[k]*sig2[k])/sig1[k]);}
        for (int k=0; k<n; k++){
            w1.push_back(0);
            w1p.coeffRef(k)=0;
            for (SparseMatrix<double,ColMajor>::InnerIterator it(A,k);it;++it){
                w1[k] += -(1/tau)*A.coeffRef(it.row(),it.col())*(-1/fu1[it.row()] + 1/fu2[it.row()]);
                w1p.coeffRef(k)+= - A.coeffRef(it.row(),it.col())*((sig2[it.row()]/sig1[it.row()])*w2[it.row()]);}
            w1p.coeffRef(k)+= w1[k];}
        
        /* H11p */
        MatXd H11p(n,n);
        SpMat AtAt = AtA.transpose();
        for (int k=0; k<n*n; k++){
            sigy.push_back(0);
            for (SparseMatrix<double,ColMajor>::InnerIterator it(AtAt,k);it;++it){
                sigy[k]+=AtAt.coeffRef(it.row(),it.col())*sigx[it.row()];}
            H11p.coeffRef(floor(k/n),k-floor(k/n)*n)=sigy[k];}
        SpMat H11ps = H11p.sparseView();
        
        /* solve */
        ConjugateGradient<SpMat > solver(H11ps);
        //SparseQR<SpMat,COLAMDOrdering<int> > solver(H11ps);
        //solver.compute(H11ps);
        VecXd dx = solver.solve(w1p);
        //VecXd dx = H11p.householderQr().solve(w1p);
        
        /*	[dx, hcond] = linsolve(H11p, w1p);
         *  if (hcond < 1e-14)
         *      disp('Matrix ill-conditioned.  Returning previous iterate.');
         *      dist('(See Section 4 of notes for more information.)');
         *      xp = x;
         *      return;
         *  end */
        VecXd Adx = A*dx;
        for (int k=0; k<m; k++){
            du.push_back((w2[k] - sig2[k]*Adx.coeffRef(k))/sig1[k]);
            dlamu1.push_back(-(lamu1[k]/fu1[k])*(Adx[k]-du[k])-lamu1[k]-(1/tau)*(1/fu1[k]));
            dlamu2.push_back( (lamu2[k]/fu2[k])*(Adx[k]+du[k])-lamu2[k]-(1/tau)*(1/fu2[k]));}
        for (int k=0; k<n; k++){
            Atdv.push_back(0);
            for (SparseMatrix<double,ColMajor>::InnerIterator it(A,k);it;++it){
                Atdv[k]+= A.coeffRef(it.row(),it.col())*(dlamu1[it.row()]-dlamu2[it.row()]);}}
        // make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
        double s=1;
        for (int k=0; k<m; k++){
            if(dlamu1[k]<0){s=std::min(s,-lamu1[k]/dlamu1[k]);};
            if(dlamu2[k]<0){s=std::min(s,-lamu2[k]/dlamu2[k]);};
            if(( Adx[k]-du[k])>0){s=MIN(s,-fu1[k]/( Adx[k]-du[k]));};
            if((-Adx[k]-du[k])>0){s=MIN(s,-fu2[k]/(-Adx[k]-du[k]));};}
        s = 0.99*s;
        
        /* backtrack */
        bool suffdec=false;
        int backiter=0;
        double resnormp;
        while(!suffdec){
            Atvp.clear();
            fu1p.clear();
            fu2p.clear();
            lamu1p.clear();
            lamu2p.clear();
            resnormp = 0;
            for (int k=0; k<n; k++){
                xp.coeffRef(k) = x.coeffRef(k) + s*dx.coeffRef(k);
                Atvp.push_back(Atv[k] + s*Atdv[k]);
                resnormp += Atvp[k]*Atvp[k];}
            for (int k=0; k<m; k++){
                up.coeffRef(k)  = u.coeffRef(k)  + s*du[k];
                Axp.coeffRef(k) = Ax.coeffRef(k) + s*Adx.coeffRef(k);
                fu1p.push_back( Axp.coeffRef(k) - y.coeffRef(k) - up.coeffRef(k));
                fu2p.push_back(-Axp.coeffRef(k) + y.coeffRef(k) - up.coeffRef(k));
                lamu1p.push_back(lamu1[k] + s*dlamu1[k]);
                lamu2p.push_back(lamu2[k] + s*dlamu2[k]);
                temp      = -lamu1p[k]*fu1p[k]-(1/tau);
                resnormp += temp*temp;
                temp      = -lamu2p[k]*fu2p[k]-(1/tau);
                resnormp += temp*temp;
                temp      = -lamu1p[k]-lamu2p[k]+1;
                resnormp += temp*temp;}
            resnormp = sqrt(resnormp);
            suffdec = (resnormp <= (1-alpha*s)*resnorm);
            s = beta*s;
            backiter = backiter + 1;
            if(backiter>32){
                std::cout<<"Stuck backtracking, returning last iterate."<<std::endl;
                std::cout<<"See Section 4 of notes for more information."<<std::endl;
                xp = x;
                break;}}
        
        /* next iteration */
        x = xp; u = up; Ax = Axp; Atv = Atvp;
        fu1 = fu1p; fu2 = fu2p; lamu1 = lamu1p; lamu2 = lamu2p;
        
        /* clear vectors */
        w1.clear(); w2.clear();
        sig1.clear(); sig2.clear(); sigx.clear(); sigy.clear();
        Atdv.clear(); du.clear(); dlamu1.clear(); dlamu2.clear();
        
        /* surrogate duality gap */
        sdg=0;
        for (int k=0; k<m; k++){
            sdg = sdg - fu1[k]*lamu1[k] - fu2[k]*lamu2[k];}
        tau = mu*2*m/sdg;
        resnorm=0;
        for (int k=0; k<n; k++){
            resnorm += Atv[k]*Atv[k];}
        for (int k=0; k<m; k++){
            temp     = -lamu1[k]*fu1[k]-(1/tau);
            resnorm += temp*temp;
            temp     = -lamu2[k]*fu2[k]-(1/tau);
            resnorm += temp*temp;
            temp     = -lamu1p[k]-lamu2p[k]+1;
            resnorm += temp*temp;}
        resnorm = sqrt(resnorm);
        
        done = ((sdg < pdtol) | (pditer >= pdmaxiter));}
    
    return x;}
