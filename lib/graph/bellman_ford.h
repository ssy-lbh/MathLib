#ifndef BELLMAN_FORD_H
#define BELLMAN_FORD_H

void bellman_ford(int n, int m, int s, int e, int *head, int *to, int *next, int *w, int *dis, int *pre){
    for (int i = 0; i < n; i++){
        dis[i] = 0x3f3f3f3f;
    }
    dis[s] = 0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            int u = to[j];
            int v = head[j];
            int t = dis[u] + w[j];
            if (t < dis[v]){
                dis[v] = t;
                pre[v] = u;
            }
        }
    }
}

#endif /* BELLMAN_FORD_H */
