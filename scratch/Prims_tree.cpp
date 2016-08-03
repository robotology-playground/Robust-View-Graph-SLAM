#include<bits/stdc++.h>
using namespace std;

int V = 5;

struct Graph
{
    int v, w;
    vector< pair<int, int>>* arr;
};

Graph* initGraph(int v)
{
    Graph* g = new Graph;

    g->v = v;
    g->w = 0;
    g->arr = new vector< pair<int, int>>[g->v];

    return g;
}

void addEdge(Graph* g, int u, int v, int w)
{
    g->arr[u].push_back(make_pair(v, w));
    g->arr[v].push_back(make_pair(u, w));
}

int minKey(int key[], int mstSet[])
{
    int mini = INT_MAX, minInd;
    for(int i=0; i<V; i++)
        if(key[i] < mini && !mstSet[i])
        {
            mini = key[i];
            minInd = i;
        }
    return minInd;
}

void printGraph(int parent[], int key[], int n)
{
	printf("Prim's MST\n");
    for(int i=1; i<n; i++)
    {
        printf("%d---%d\t\t%d\n", parent[i], i, key[i]);
    }
}

void Prims(Graph* g)
{
    int parent[g->v];
    int mstSet[g->v];
    int key[g->v];

    for(int i=0; i< g->v; i++)
    {
        key[i] = INT_MAX;
        mstSet[i] = 0;
    }

    parent[0] = -1;
    key[0] = 0;
    
    for(int i=0; i< g->v-1; i++)
    {
        int u = minKey(key, mstSet);

        mstSet[u] = 1;

        vector< pair<int, int>>::iterator it;
        for(it=g->arr[u].begin(); it!=g->arr[u].end(); it++)
        {
            int v = (*it).first;
            int w = (*it).second;
            if(!mstSet[v] && key[v] > w)
            {
                parent[v] = u;
                key[v] = w;
            }
        }
    }

    printGraph(parent, key, g->v);
}

void printGraph(Graph* g)
{
	printf("Original Graph as Adj. List\n");
    for(int i=0; i<g->v; i++)
    {
        printf("%d-> ", i);

        vector< pair<int, int>>::iterator it;
        for(it=g->arr[i].begin(); it!=g->arr[i].end(); it++)
        {
            printf("(%d, %d) ", it->first, it->second);
        }
        printf("\n");
    }
}

int main()
{
    /** Let us create the following graph
          2    3
      (0)--(1)--(2)
       |   / \   |
      6| 8/   \5 |7
       | /     \ |
      (3)-------(4)
            9          */
    Graph* g = initGraph(V);

    addEdge(g, 0, 1, 2);
    addEdge(g, 1, 2, 3);
    addEdge(g, 2, 4, 7);
    addEdge(g, 4, 1, 5);
    addEdge(g, 4, 3, 9);
    addEdge(g, 3, 1, 8);
    addEdge(g, 3, 0, 6);
    printGraph(g);printf("\n");

    Prims(g);

    return 0;
}

