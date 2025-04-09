#ifndef HYPERIBLT_H
#define HYPERIBLT_H

#include <bits/stdc++.h>
#include "tools.cpp"
using namespace std;

struct Cell {
    int c; poly p; Cell();
} ;
struct HyperIBLT {
    vector<Cell> B;
    int size();
    void init();
    Cell & operator [] (int);
    friend HyperIBLT operator - (HyperIBLT, HyperIBLT);
    void InsertToCell(int, int);
    void Update(int);
    pair<vector<int>, vector<int> > PureCellDecode(int);
    bool PureCellVerify(int);
    void Extract(int, int);
    variant<pair<vi, vi>, bool> Decode();
} ;
HyperIBLT Encode(vector<int>) ;
#endif