#include <bits/stdc++.h>
#include "HyperIBLT.cpp"

using namespace std;
#define fo(v, a, b) for(int v = a; v <= b; v++)
#define fr(v, a, b) for(int v = a; v >= b; v--)
#define cl(a, v) memset(a, v, sizeof(a))

const int MAXLEN = (1 << 19);

vector<int> AliceData, BobData; vector<int> Adiff, Bdiff;

int D;
void genereteData(int cA, int cB) {
    assert(D >= abs(cA - cB));
    AliceData.resize(max(cA, cB));
    set<int> S, Sp;
    for(int i = 0; i < max(cA, cB); i++) {
        int x = rng() % (P - 1) + 1;
        while(S.find(x) != S.end()) x = rng() % (P - 1) + 1;
        AliceData[i] = x, S.insert(x);
    }
    BobData = AliceData;
    AliceData.resize(cA), BobData.resize(cB);
    for(int i = 1; i <= (D - abs(cA - cB)) / 2; i++) {
        int p = rng() % min(cA, cB), v = rng() % (P - 1) + 1;
        while(Sp.find(p) != Sp.end()) p = rng() % min(cA, cB);
        while(S.find(v) != S.end()) v = rng() % (P - 1) + 1;
        BobData[p] = v, S.insert(v), Sp.insert(p);
        Adiff.push_back(AliceData[p]), Bdiff.push_back(BobData[p]);
    }
    for(int i = cB; i < cA; i++) Adiff.push_back(AliceData[i]);
    for(int i = cA; i < cB; i++) Bdiff.push_back(BobData[i]);
    sort(Adiff.begin(), Adiff.end()), sort(Bdiff.begin(), Bdiff.end());
    shuffle(AliceData.begin(), AliceData.end(), rng);
    shuffle(BobData.begin(), BobData.end(), rng);
}

int main() {
    tool::init(MAXLEN), tool::Pinit(MAXLEN);

    int cA, cB;
    // cA = 10, cB = 10, D = 13;
    // cA = 1000000, cB = 1000000, D = 10000;
    // cA = 3000, cB = 3000, D = 2000;
    // cA = 1000000, cB = 1000000, D = 1000;
    // cA = 10000000, cB = 10000000, D = 10000;
    cA = 500000, cB = 500000, D = 100000;
    // scanf("%d%d%d", &cA, &cB, &D);
    genereteData(cA, cB);
    // AliceData = {233613855, 84053574, 955552696, 443195438};
    // BobData = {955552696, 233613855, 84053574, 356042544, 443195438};
    // print(AliceData), print(BobData);
    // print(Adiff), print(Bdiff);

    int Znow;
    // l = 5, k = 2, M = 20, Znow = 0;
    // l = 5, k = 2, M = 2400, Znow = 4;
    // l = 5, k = 2, M = 1200, Znow = 2;
    // l = 6, k = 2, M = 11500, Znow = 5;
    l = 6, k = 2, M = 107500, Znow = 10;
    // M = 108500;
    M /= l, HashingInit(Znow);
    // printf("l = %d k = %d M = %d\n", l, k, M);
    // for(int x : AliceData) fo(j, 1, k)
    //     printf("+ %d %d : %d\n", j, x, h(j, x));
    // for(int x : BobData) fo(j, 1, k)
    //     printf("- %d %d : %d\n", j, x, h(j, x));
    auto check = [&](Cell &v) {
        return 0 <= v.c && v.c <= 2 * l && v.p.size() <= l;
    } ;
    double EncodeBegin = clock();
    cerr << "BeginAliceEncoding" << endl;
        auto Alice = Encode(AliceData);
    double EncodeEnd = clock();
    cerr << "BeginBobEncoding" << endl;
        auto Bob = Encode(BobData);
    double DecodeBegin = clock();
    cerr << "BeginDecoding" << endl;
        auto DiffRes = (Alice - Bob).Decode();
    double DecodeEnd = clock();

    if(DiffRes.index() == 1) puts("ERROR"), exit(0);
    auto Diff = get<pair<vector<int>, vector<int> > >(DiffRes);
    auto ansad = Diff.first, ansbd = Diff.second;
    assert(ansad.size() == Adiff.size() && ansbd.size() == Bdiff.size());
    for(int i = 0; i < (int)Adiff.size(); i++) assert(Adiff[i] == ansad[i]);
    for(int i = 0; i < (int)Bdiff.size(); i++) assert(Bdiff[i] == ansbd[i]);
    puts("AC");
    // puts("Alice - Bob : ");
    // for(int x : Diff.first) printf("%d ", x);
    // puts("");
    // puts("Bob - Alice : ");
    // for(int x : Diff.second) printf("%d ", x);
    // puts("");

    printf("Running Time of Encode : %.3lf s\n", (EncodeEnd - EncodeBegin) / CLOCKS_PER_SEC);
    printf("Running Time of Decode : %.3lf s\n", (DecodeEnd - DecodeBegin) / CLOCKS_PER_SEC);
    printf("Com Cost : %d\n", Bob.queryMemory());

    return 0;
}
/*
Sample :
Input:
10 3
1 2 3 4 5 6 7 8 9 10
1 12 3 4 5 6 7 28 39 10
Output:
Alice - Bob :
2 8 9
Bob - Alice :
12 28 39
*/