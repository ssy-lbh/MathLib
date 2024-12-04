#ifndef STRUCTURE_TREE_ARRAY_H
#define STRUCTURE_TREE_ARRAY_H

#include "math/tensor.h"

// 高维树状数组，单点修改，区间查询

template <typename T, int L>
void tree_array_add(TTensor<T, L>& tree, int i, T x){
    for (; i < L; i += i & -i)
        tree[i] += x;
}

template <typename T, int IL>
void tree_array_add(T& tree, TTensor<int, IL>& idx, T x, int depth = 0){
    tree += x;
}

template <typename T, int CL, int... L, int IL>
void tree_array_add(TTensor<T, CL, L...>& tree, TTensor<int, IL>& idx, T x, int depth = 0){
    static_assert(sizeof...(L) < IL, "tree_array_add: index length mismatch");
    for (int i = idx[depth]; i < CL; i += i & -i)
        tree_array_add(tree[i], idx, x, depth + 1);
}

template <typename T, int L>
T tree_array_sum(TTensor<T, L>& tree, int i){
    T res = zero(T());
    for (; i; i -= i & -i)
        res += tree[i];
    return res;
}

template <typename T, int IL>
T tree_array_sum(T& tree, TTensor<int, IL>& idx, int depth = 0){
    return tree;
}

template <typename T, int CL, int... L, int IL>
T tree_array_sum(TTensor<T, CL, L...>& tree, TTensor<int, IL>& idx, int depth = 0){
    static_assert(sizeof...(L) < IL, "tree_array_sum: index length mismatch");
    T res = zero(T());
    for (int i = idx[depth]; i; i -= i & -i)
        res += tree_array_sum(tree[i], idx, depth + 1);
    return res;
}

#endif
