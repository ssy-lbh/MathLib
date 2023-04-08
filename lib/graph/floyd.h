#ifndef FLOYD_H
#define FLOYD_H

template <typename T, int N>
void floyd(T (&a)[N][N]){
    for (int k = 0; k < N; k++){
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                T t = a[i][k] + a[k][j];
                a[i][j] = t < a[i][j] ? t : a[i][j];
            }
        }
    }
}

#endif /* FLOYD_H */
