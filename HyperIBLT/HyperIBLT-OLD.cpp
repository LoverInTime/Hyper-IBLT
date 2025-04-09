#include <bits/stdc++.h>
#include "tools.cpp"
#include "hash.cpp"
using namespace std;
#define fo(v, a, b) for(int v = a; v <= b; v++)
#define fr(v, a, b) for(int v = a; v >= b; v--)
#define cl(a, v) memset(a, v, sizeof(a))

// using namespace RandomHash;
using namespace SpatialCoupling;

int k, l, d;
queue<int> Q; vector<bool> Vis;
struct Cell {
    short c; poly p;
    Cell() {  c = 0, p = plv({1});  }
} ;
struct HyperIBLT {
    vector<Cell> B;
    int queryMemory() {  return B.size() * (sizeof(short) + B[0].p.size() * sizeof(int));  }
    int size() {  return B.size();  }
    void init() {  B.resize(M); fo(i, 0, M - 1) B[i].p.rs(l);  }
    Cell & operator [] (int x) { return B[x]; }
    friend HyperIBLT operator - (HyperIBLT a, HyperIBLT b) {
        HyperIBLT res; res.init();
        fo(i, 0, M - 1) {
            res[i].c = (a[i].c - b[i].c + (2 * l + 1)) % (2 * l + 1);
            a[i].p.rs(l), b[i].p.rs(l);
            res[i].p = a[i].p * b[i].p.Inv(), res[i].p.rs(l);
        }
        return res;
    }

    void InsertToCell(int i, int x) {
        B[i].c = (B[i].c + 1) % (2 * l + 1);
        fr(j, l - 1, 1) B[i].p[j] = ((ll)(P - x) * B[i].p[j] + B[i].p[j - 1]) % P;
        B[i].p[0] = (ll)(P - x) * B[i].p[0] % P;
        // B[i].p = B[i].p * plv({P - x, 1}), B[i].p.rs(l);
    }
    void Update(int x) {
        fo(i, 1, k) InsertToCell(h(i, x), x);
    }

    pair<vector<int>, vector<int> > PureCellDecode(int i) {
        int m = B[i].c; if(m > l) m -= 2 * l + 1;
        auto ChiRes = tool :: RFuncReconstruct(B[i].p, l, m);
        if(ChiRes.index() == 1) assert(0);
        auto Chi = get<pair<poly, poly> >(ChiRes);
        auto DeltaA = tool :: findRoots(Chi.first), DeltaB = tool :: findRoots(Chi.second);
        sort(DeltaA.begin(), DeltaA.end()), sort(DeltaB.begin(), DeltaB.end());
        return make_pair(DeltaA, DeltaB);
    }
    bool PureCellVerify(int i) {
        int m = B[i].c; if(m > l) m -= 2 * l + 1;
        auto ChiRes = tool :: RFuncReconstruct(B[i].p, l, m);
        if(ChiRes.index() == 1) return false;
        auto Chi = get<pair<poly, poly> >(ChiRes);
        auto DA = tool :: findAllRoots(Chi.first), DB = tool :: findAllRoots(Chi.second);
        if(DA.index() == 1 || DB.index() == 1) return false;
        auto DeltaA = get<vi>(DA), DeltaB = get<vi>(DB);
        for(auto t : {DeltaA, DeltaB}) for(int x : t) {
            bool fl = false;
            fo(j, 1, k) if(h(j, x) == i) {  fl = true; break;  }
            if(!fl) return false;
        }
        return true;
    }
    void Extract(int x, int type) {
        poly now = plv({P - x, 1});
        if(type == 0) now.rs(l), now = now.Inv();
        fo(j, 1, k) {
            int i = h(j, x);
            if(type == 0)
                B[i].c = (B[i].c + 2 * l) % (2 * l + 1);
            else
                B[i].c = (B[i].c + 1) % (2 * l + 1);
            B[i].p = B[i].p * now, B[i].p.rs(l);
            if(!Vis[i] && PureCellVerify(i))
                Vis[i] = true, Q.push(i);
        }
    }

    variant<pair<vi, vi>, bool> Decode() {
        while(Q.size()) Q.pop();
        Vis.clear(), Vis.resize(M, false);
        vector<int> SA, SB;
        fo(i, 0, M - 1) if(PureCellVerify(i))
            Vis[i] = true, Q.push(i);
        while(Q.size()) {
            int i = Q.front(); Q.pop();
            auto [da, db] = PureCellDecode(i);
            if(da.empty() && db.empty()) continue;
            for(int x : da) Extract(x, 0), SA.push_back(x);
            for(int x : db) Extract(x, 1), SB.push_back(x);
        }
        fo(i, 0, M - 1) if(!Vis[i])
            return false;
        sort(SA.begin(), SA.end()), sort(SB.begin(), SB.end());
        return make_pair(SA, SB);
    }
} ;

HyperIBLT Encode(vector<int> v) {
    HyperIBLT res; res.init();
    for(int x : v) res.Update(x);
    return res;
}