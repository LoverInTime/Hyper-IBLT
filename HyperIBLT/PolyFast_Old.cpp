#include <bits/stdc++.h>
using namespace std;
#define fo(v, a, b) for(int v = a; v <= b; v++)
#define fr(v, a, b) for(int v = a; v >= b; v--)
#define cl(a, v) memset(a, v, sizeof(a))

const int P = 998244353;
typedef vector<int> vi;
mt19937 rng(time(0));

namespace tool {
    // #define fo(v, a, b) for(int v = a; v <= b; v++)
    // #define fr(v, a, b) for(int v = a; v >= b; v--)
    // #define cl(a, v) memset(a, v, sizeof(a))
    // const int P = 998244353;
    //////////
    const int mod = 998244353, _G = 3, N = (1 << 21), inv2 = (mod + 1) / 2;
    #define sz(a) ((int)a.size())
    #define L(i, j, k) for(int i = (j); i <= (k); i++)
    #define R(i, j, k) for(int i = (j); i >= (k); i--)
    #define ll long long
    #define vi vector<int>
    #define add(a, b) (a + b >= mod ? a + b - mod : a + b)
    #define dec(a, b) (a < b ? a - b + mod : a - b)
    int qpow(int x, int y = mod - 2) {
        int res = 1;
        for(; y; x = (ll) x * x % mod, y >>= 1) if(y & 1) res = (ll) res * x % mod;
        return res;
    }
    int fac[N + 1], ifac[N + 1], inv[N + 1];
    void init(int x) {
        fac[0] = ifac[0] = inv[1] = 1;
        L(i, 2, x) inv[i] = (ll) inv[mod % i] * (mod - mod / i) % mod;
        L(i, 1, x) fac[i] = (ll) fac[i - 1] * i % mod, ifac[i] = (ll) ifac[i - 1] * inv[i] % mod;
    }
    int C(int x, int y) {
        return y < 0 || x < y ? 0 : (ll) fac[x] * ifac[y] % mod * ifac[x - y] % mod;
    }
    inline int sgn(int x) {
        return (x & 1) ? mod - 1 : 1;
    }
    int rt[N], Lim;
    void Pinit(int x) {
        for(Lim = 1; Lim <= x; Lim <<= 1) ;
        for(int i = 1; i < Lim; i <<= 1) {
            int sG = qpow (_G, (mod - 1) / (i << 1));
            rt[i] = 1;
            L(j, i + 1, i * 2 - 1) rt[j] = (ll) rt[j - 1] * sG % mod;
        }
    }
    struct poly {
        vector<int> a;
        void print() {
            for(int i = 0; i < sz(a); i++) printf("%d ", a[i]);
            puts("");
        }
        void read(int n) {
            a.resize(n);
            for(int i = 0; i < n; i++) scanf("%d", &a[i]);
        }
        int size() { return sz(a); }
        void PopZero() {
            while(a.size() && a.back() == 0) a.pop_back();
        }
        int deg() {
            PopZero(); return sz(a) == 1 && a[0] == 0 ? -1e9 : sz(a) - 1;
        }
        int & operator [] (int x) { return a[x]; }
        int v(int x) { return x < 0 || x >= sz(a) ? 0 : a[x]; }
        void clear() { vector<int> ().swap(a); }
        void rs(int x = 0) { a.resize(x); }
        poly (unsigned int n = 0) { rs(n); }
        poly (vector<int> o) { a = o; }
        // poly (const poly &o) { a = o.a; }
        poly Rs(int x = 0) { vi res = a; res.resize(x); return res; }
        inline void dif() {
            int n = sz(a);
            for (int l = n >> 1; l >= 1; l >>= 1) 
                for(int j = 0; j < n; j += l << 1) 
                    for(int k = 0, *w = rt + l; k < l; k++, w++) {
                        int x = a[j + k], y = a[j + k + l];
                        a[j + k] = add(x, y);
                        a[j + k + l] = (ll) * w * dec(x, y) % mod;
                    }
        }
        void dit () {
            int n = sz(a);
            for(int i = 2; i <= n; i <<= 1) 
                for(int j = 0, l = (i >> 1); j < n; j += i) 
                    for(int k = 0, *w = rt + l; k < l; k++, w++) {
                        int pa = a[j + k], pb = (ll) a[j + k + l] * *w % mod;
                        a[j + k] = add(pa, pb), a[j + k + l] = dec(pa, pb);
                    }
            reverse(a.begin() + 1, a.end());
            for(int i = 0, iv = qpow(n); i < n; i++) a[i] = (ll) a[i] * iv % mod;
        } 
        friend poly operator * (poly aa, poly bb) {
            if(!sz(aa) || !sz(bb)) return (vi){};
            int lim, all = sz(aa) + sz(bb) - 1;
            for(lim = 1; lim < all; lim <<= 1);
            aa.rs(lim), bb.rs(lim), aa.dif(), bb.dif();
            L(i, 0, lim - 1) aa[i] = (ll) aa[i] * bb[i] % mod;
            aa.dit(), aa.a.resize(all);
            return aa;
        }
        poly Inv() {
            poly res, f, g;
            res.rs(1), res[0] = qpow(a[0]);
            for(int m = 1, pn; m < sz(a); m <<= 1) {
                pn = m << 1, f = res, g.rs(pn), f.rs(pn);
                for(int i = 0; i < pn; i++) g[i] = (*this).v(i);
                f.dif(), g.dif();
                for(int i = 0; i < pn; i++) g[i] = (ll) f[i] * g[i] % mod;
                g.dit();
                for(int i = 0; i < m; i++) g[i] = 0;
                g.dif();
                for(int i = 0; i < pn; i++) g[i] = (ll) f[i] * g[i] % mod;
                g.dit(), res.rs(pn);
                for(int i = m; i < min(pn, sz(a)); i++) res[i] = (mod - g[i]) % mod;
            } 
            return res.rs(sz(a)), res;
        }
        poly Shift (int x) {
            assert(sz(a) + x > 0);
            poly zm (sz(a) + x);
            L(i, max(-x, 0), sz(a) - 1) zm[i + x] = a[i];
            return zm; 
        }
        void ShiftSelf (int x) {
            assert(sz(a) + x > 0); vi zm (sz(a) + x);
            L(i, max(-x, 0), sz(a) - 1) zm[i + x] = a[i];
            a = zm;
        }
        friend poly operator * (poly aa, int bb) {
            poly res(sz(aa));
            L(i, 0, sz(aa) - 1) res[i] = (ll) aa[i] * bb % mod;
            return res;
        }
        friend poly operator * (int bb, poly aa) {
            return aa * bb;
        }
        friend bool operator == (poly aa, poly bb) {
            aa.PopZero(), bb.PopZero();
            if(aa.size() != bb.size()) return false;
            L(i, 0, sz(aa) - 1) if(aa[i] != bb[i]) return false;
            return true;
        }
        friend poly operator + (poly aa, poly bb) {
            vector<int> res(max(sz(aa), sz(bb)));
            L(i, 0, sz(res) - 1) res[i] = add(aa.v(i), bb.v(i));
            return poly(res);
        }
        friend poly operator - (poly aa, poly bb) {
            vector<int> res(max(sz(aa), sz(bb)));
            L(i, 0, sz(res) - 1) res[i] = dec(aa.v(i), bb.v(i));
            return poly(res);
        }
        poly & operator += (poly o) {
            rs(max(sz(a), sz(o)));
            L(i, 0, sz(a) - 1) (a[i] += o.v(i)) %= mod;
            return (*this);
        }
        poly & operator -= (poly o) {
            rs(max(sz(a), sz(o)));
            L(i, 0, sz(a) - 1) (a[i] += mod - o.v(i)) %= mod;
            return (*this);
        }
        poly & operator *= (poly o) {
            return (*this) = (*this) * o;
        }
        poly Integ() {
            if(!sz(a)) return poly();
            poly res(sz(a) + 1);
            L(i, 1, sz(a)) res[i] = (ll) a[i - 1] * inv[i] % mod;
            return res;
        }
        poly Deriv() {
            if(!sz(a)) return poly();
            poly res(sz(a) - 1); 
            L(i, 1, sz(a) - 1) res[i - 1] = (ll) a[i] * i % mod;
            return res;
        }
        poly Ln() {
            poly g = ((*this).Inv() * (*this).Deriv()).Integ();
            return g.rs(sz(a)), g;
        }
        poly Exp() {
            poly res(1), f; 
            res[0] = 1;
            for(int m = 1, pn; m < sz(a); m <<= 1) {
                pn = min(m << 1, sz(a)), f.rs(pn), res.rs(pn);
                for(int i = 0; i < pn; i++) f[i] = (*this).v(i);
                f -= res.Ln(), (f[0] += 1) %= mod, res *= f, res.rs(pn); 
            }
            return res.rs(sz(a)), res;
        }
        poly pow(int x, int rx = -1) { // x : the power % mod; rx : the power % (mod - 1)
            if(rx == -1) rx = x;
            int cnt = 0;
            while (a[cnt] == 0 && cnt < sz(a)) cnt += 1;
            
            poly res = (*this);
            L(i, cnt, sz(a) - 1) res[i - cnt] = res[i];
            L(i, sz(a) - cnt, sz(a) - 1) res[i] = 0;
            int c = res[0], w = qpow (res[0]);
            L(i, 0, sz(res) - 1) res[i] = (ll) res[i] * w % mod;
            res = res.Ln();
            L(i, 0, sz(res) - 1) res[i] = (ll) res[i] * x % mod;
            res = res.Exp();
            c = qpow (c, rx);
            L(i, 0, sz(res) - 1) res[i] = (ll) res[i] * c % mod;
            
            if((ll) cnt * x > sz(a)) L(i, 0, sz(a) - 1) res[i] = 0;
            else if(cnt) {
                R(i, sz(a) - cnt * x - 1, 0) res[i + cnt * x] = res[i];
                L(i, 0, cnt * x - 1) res[i] = 0; 
            }
            return res;
        }
        poly sqrt(int Rt = 1) {
            poly res(1), f; 
            res[0] = Rt;
            for(int m = 1, pn; m < sz(a); m <<= 1) {
                pn = min(m << 1, sz(a)), f.rs(pn);
                for(int i = 0; i < pn; i++) f[i] = (*this).v(i);
                f += res * res, f.rs(pn), res.rs(pn), res = f * res.Inv(), res.rs(pn);
                for(int i = 0; i < pn; i++) res[i] = (ll) res[i] * inv2 % mod;
            } 
            return res;
        }
        friend poly mul (poly aa, poly bb, int k) {
            if(!sz(aa) || !sz(bb)) return {};
            int lim; 
            for(lim = 1; lim < k; lim <<= 1);
            aa.rs(lim), bb.rs(lim), aa.dif(), bb.dif();
            L(i, 0, lim - 1) aa[i] = (ll) aa[i] * bb[i] % mod;
            aa.dit(), aa.a.resize(lim);
            return aa;
        }
        void Rev() {
            reverse(a.begin(), a.end());
        }
        friend pair < poly, poly > div (poly f, poly g) { /* f / g = first, f % g = second */
            f.rs(max(sz(f), sz(g))), f.Rev(), g.Rev();
            int n = sz(f), m = sz(g);
            poly A = g.Rs(n - m + 1).Inv(), t;
            A *= f.Rs(n - m + 1), A.rs(n - m + 1), A.Rev(), g.Rev(), f.Rev(), t = f - A * g, t.rs(m - 1);
            A.PopZero(), t.PopZero(); return make_pair(A, t);
        }
        void monic() {
            assert(a.size()); auto C = qpow(a.back());
            for(int &x : a) x = (ll)x * C % P;
        }
    } ;
    inline poly plv(vi v) {  return poly(v);  }
    struct polyMat {
        poly a00, a01, a10, a11;
        polyMat operator*(const polyMat &t) const {
            polyMat res;
            res.a00 = a00 * t.a00 + a01 * t.a10, res.a01 = a00 * t.a01 + a01 * t.a11;
            res.a10 = a10 * t.a00 + a11 * t.a10, res.a11 = a10 * t.a01 + a11 * t.a11;
            res.a00.PopZero(), res.a01.PopZero(), res.a10.PopZero(), res.a11.PopZero();
            return res;
        }
    } ;
    inline polyMat gen(poly Q) {
        return {plv({0}), plv({1}), plv({1}), (P - 1) * Q};
    }
    inline void Mul(polyMat M, poly &A, poly &B) {
        poly tA = A, tB = B;
        A = M.a00 * tA + M.a01 * tB, B = M.a10 * tA + M.a11 * tB;
        A.PopZero(), B.PopZero(); 
    }
    inline pair<poly, poly> Mul(polyMat M, pair<poly, poly> u) {
        Mul(M, u.first, u.second); return u;
    }
    polyMat HalfGCD(poly A, poly B) {
        if(B.deg() == 0) return {plv({1}), plv({0}), plv({0}), plv({1})};
        if(A.deg() == 0) return {plv({0}), plv({1}), plv({1}), plv({0})};
        int len = A.deg(), m = (len + 1) / 2;
        if(B.deg() < m) return {plv({1}), plv({0}), plv({0}), plv({1})};
        poly A1 = A.Shift(-m), B1 = B.Shift(-m); auto M = HalfGCD(A1, B1);
        Mul(M, A, B); if(B.deg() < m) return M;

        auto tmp = div(A, B); A = B, B = tmp.second;
        M = gen(tmp.first) * M; if(B.deg() < m) return M;

        int k = 2 * m - A.deg(); A1 = A.Shift(-k), B1 = B.Shift(-k);
        return HalfGCD(A1, B1) * M;
    }
    polyMat coGCD(poly A, poly B) {
        auto M = HalfGCD(A, B); Mul(M, A, B);
        if(B.size() == 0) return M;
        auto tmp = div(A, B); A = B, B = tmp.second;
        M = gen(tmp.first) * M; if(B.size() == 0) return M;
        return coGCD(A, B) * M;
    }
    inline poly GCD(poly A, poly B) {
        if(B == plv({0})) return A;
        if(A == plv({0})) return B;
        auto M = coGCD(A, B); auto res = M.a00 * A + M.a01 * B;
        res.PopZero(), res.monic();
        return res;
    }
    // deg A < 2 * n + 1, deg p, deg q <= n
    pair<poly, poly> Recon(poly A, int n) {
        poly M(2 * n + 2); M[2 * n + 1] = 1; A.rs(2 * n + 1);
        auto Mat = HalfGCD(A, M);
        auto X = Mat.a00 * A + Mat.a01 * M;
        X.PopZero(), Mat.a00.PopZero();
        if(X.deg() <= n && Mat.a00.deg() <= n)
            return make_pair(X, Mat.a00);
        auto Y = Mat.a10 * A + Mat.a11 * M;
        Y.PopZero(), Mat.a10.PopZero();
        if(Y.deg() <= n && Mat.a10.deg() <= n)
            return make_pair(Y, Mat.a10);
        assert(0);
    }
    // find (p, q) s.t.
    // p / q mod x^n = A, deg p + deg q < n, deg p - deg q = m, [x^0]p = [x^0]q = 1
    // false for no solution
    variant<pair<poly, poly>, bool> Reconstruct(poly A, int n, int m) {
        int pdeg = (n - 1 + m) / 2, qdeg = (n - 1 - m) / 2;
        assert(pdeg <= qdeg); int delta = qdeg - pdeg; A.ShiftSelf(delta);
        auto res = Recon(A, qdeg);
        bool fl = (res.first.deg() >= delta);
        for(int i = 0; fl && i < delta; i++)
            fl = (res.first[i] == 0);
        if(!fl) return false;
        return make_pair(res.first.Shift(-delta), res.second);
    }
    variant<pair<poly, poly>, bool> RFuncReconstruct(poly A, int n, int m) {
        A.rs(n); variant<pair<poly, poly>, bool> Res;
        if(m <= 0) Res = Reconstruct(A, n, m);
        else Res = Reconstruct(A.Inv(), n, -m);
        if(Res.index() == 1) return false;
        auto res = get<pair<poly, poly> >(Res);
        if(m > 0) swap(res.first, res.second);
        auto fz = res.first.Rs(n), fm = res.second.Rs(n);
        auto now = fz * fm.Inv(); now.rs(n);
        if(now == A) return res;
        return false;
    }

    inline poly randPoly(int n) { // deg = n - 1
        poly r(n);
        for(int i = 0; i < n; i++) r[i] = rng() % P;
        return r;
    }
    poly power(poly A, int n, poly M) {
        auto res = plv({1}), a = A;
        while(n) {
            if(n & 1) res = div(res * a, M).second;
            a = div(a * a, M).second, n >>= 1;
        }
        return res;
    }
    // vector<int> divide(poly A) {
    //     if(A.deg() == 0) return {};
    //     if(A.deg() == 1) return {P - A[0]};
    //     poly R = power(randPoly(A.size()), (P - 1) >> 1, A);
    //     auto lp = GCD(R - plv({1}), A), rp = div(A, lp).first;
    //     auto lv = divide(lp), rv = divide(rp);
    //     for(int x : rv) lv.push_back(x);
    //     return lv;
    // }
    vector<int> divide2(poly iA) {
        queue<poly> Q; Q.push(iA); vector<int> vec;
        int cnt = 0;
        while(Q.size()) {
            poly A = Q.front(); Q.pop(); cnt++;
            // printf("%d\n", A.size());
            if(A.deg() == 0) continue;
            if(A.deg() == 1) {  vec.push_back(P - A[0]); continue;  }
            poly R = power(randPoly(A.size()), (P - 1) >> 1, A);
            auto lp = GCD(R - plv({1}), A), rp = div(A, lp).first;
            Q.push(lp), Q.push(rp);
        }
        // printf("cnt = %d\n", cnt);
        return vec;
    }
    vector<int> findRoots(poly A) {
        A.monic();
        // the following codes aims to calculate gcd(X^q - X, A),
        // but since repeated roots will not occur in A, these codes is not necessary to use.
        auto UniqueRoot = [&](int q) {
            auto res = power(plv({0, 1}), q, A);
            res = res - plv({0, 1}), A = GCD(A, res);
        } ;
        UniqueRoot(P);
        auto res = divide2(A);
        sort(res.begin(), res.end());
        return res;
    }

    namespace eval {
        poly A[N], B[N], a, sav;
        int X[N], Y[N];
        void Divide1(int id, int l, int r) {
            if(l == r) return A[id] = poly(vi{1, (mod - X[l]) % mod}), void();
            int mid = (l + r) >> 1;
            Divide1(id << 1, l, mid), Divide1(id << 1 | 1, mid + 1, r);
            A[id] = A[id << 1] * A[id << 1 | 1];
        }
        void Divide2(int id, int l, int r) {
            if(l == r) return Y[l] = (a[0] + (ll) X[l] * B[id][0] % mod) % mod, void();
            int mid = (l + r) >> 1, len = r - l + 1, la = mid - l + 1, lb = r - mid;
            sav = mul(B[id], A[id << 1 | 1], len), B[id << 1].rs(la);
            L(i, 0, la - 1) B[id << 1][i] = sav[i + len - la];
            sav = mul(B[id], A[id << 1], len), B[id << 1 | 1].rs(lb);
            L(i, 0, lb - 1) B[id << 1 | 1][i] = sav[i + len - lb];
            Divide2(id << 1, l, mid), Divide2(id << 1 | 1, mid + 1, r);
        }
        vector<int> solve (poly F, vector<int> x) {
            a = F; int m = x.size();
            L(i, 1, m) X[i] = x[i - 1];
            if(sz(a) < m + 1) a.rs(m + 1);
            int n = sz(a);
            Divide1(1, 1, m);
            sav = a, sav.Rev(), sav *= A[1].Rs(n).Inv(), B[1].rs(m);
            L(i, 0, m - 1) B[1][i] = sav[i + n - m - 1];
            Divide2(1, 1, m);
            vector<int> y(m); L(i, 1, m) y[i - 1] = Y[i];
            return y;
        }
    }
}
using namespace tool;
// using tool::poly; using tool::plv;

void print(vector<int> v) {
    for(int x : v) printf("%d ", x);
    puts("");
}

poly dataToPoly(vector<int> &data, int l, int r, int Size) { // mod x^{Size}
    if(l > r) return plv({1});
    if(l == r) return plv({1, data[l]});
    int mid = (l + r) >> 1;
    poly res = dataToPoly(data, l, mid, Size) * dataToPoly(data, mid + 1, r, Size);
    return res.Rs(min(res.size(), Size));
}
inline poly CalcCharPoly(vector<int> &data, int Size) {
    return dataToPoly(data, 0, data.size() - 1, Size);
}

const int MAXLEN = (1 << 19);
// const int D = 6;
// const vector<int> AliceData = {1, 2, 3, 4}, BobData = {2, 6, 7, 5};
// const int D = 6;
// const vector<int> AliceData = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
//                   BobData = {1, 12, 3, 4, 5, 6, 7, 28, 10};
int cA, cB, D; vector<int> AliceData, BobData;

// D denotes the number of elements with different occurences in A and B
// e.g. A = {1, 2, 3}, B = {2, 3, 3} gives D = 2, since |{1, 3}| = 2 

struct info {
    int sz; poly p;
    info() {  sz = 0, p = plv({1});  }
} ;
vector<int> S;

#define len (D + 1)
void Update(int x, info &Now, vector<int> &S) {
    assert(x != 0);
    Now.sz = (Now.sz + 1) % (2 * D + 1), S.push_back(x);
    if(S.size() == D) {
        Now.p = Now.p * CalcCharPoly(S, len);
        Now.p.rs(len), S = {};
    }
}
info Encode(vector<int> Data) {
    info res; vector<int> S;
    for(int x : Data) Update(x, res, S);
    if(!S.empty())
        res.p = res.p * CalcCharPoly(S, len), S = {};
    res.p.rs(len), res.p.ShiftSelf(-1); return res;
}
variant<pair<vector<int>, vector<int> >, bool> Decode(info Alice, info Bob) {
    Alice.p.ShiftSelf(1), Alice.p[0] = 1, Bob.p.ShiftSelf(1), Bob.p[0] = 1;
    poly R = Alice.p * Bob.p.Inv(); R.rs(len);
    int m = (Alice.sz - Bob.sz + 2 * D + 1) % (2 * D + 1);
    if(m > D) m -= 2 * D + 1;
    auto ChiRes = tool :: RFuncReconstruct(R, D + 1, m);
    if(ChiRes.index() == 1) return false;
    auto Chi = get<pair<poly, poly> >(ChiRes);
    auto DeltaA = tool :: findRoots(Chi.first),
         DeltaB = tool :: findRoots(Chi.second);
    if(DeltaA.size() < Chi.first.deg() || DeltaB.size() < Chi.second.deg())
        return false;
    auto solve = [&](vi v) {
        vector<int> ans;
        for(int x : v) ans.push_back((P - qpow(x) % P) % P);
        sort(ans.begin(), ans.end());
        return ans;
    };
    return make_pair(solve(DeltaA), solve(DeltaB));
}

int main() {
    tool::init(MAXLEN), tool::Pinit(MAXLEN);

    // int type = 1; scanf("%d", &type);
    scanf("%d%d%d", &cA, &cB, &D);
    assert(D >= abs(cA - cB));
    // if((D - abs(cA - cB)) & 1) D--;
    // vector<int> AliceData(count), BobData(count);
    // if(type == 1) {
    //     for(int i = 0; i < count; i++) scanf("%d", &AliceData[i]);
    //     for(int i = 0; i < count; i++) scanf("%d", &BobData[i]);
    // } else {
        // generate data
        AliceData.resize(max(cA, cB));
        set<int> S, Sp;
        for(int i = 0; i < max(cA, cB); i++) {
            int x = rng() % (P - 1) + 1;
            while(S.find(x) != S.end()) x = rng() % (P - 1) + 1;
            AliceData[i] = x, S.insert(x);
        }
        BobData = AliceData;
        AliceData.resize(cA), BobData.resize(cB);
        vector<int> ad, bd;
        for(int i = 1; i <= (D - abs(cA - cB)) / 2; i++) {
            int p = rng() % min(cA, cB), v = rng() % (P - 1) + 1;
            while(Sp.find(p) != Sp.end()) p = rng() % min(cA, cB);
            while(S.find(v) != S.end()) v = rng() % (P - 1) + 1;
            BobData[p] = v, S.insert(v), Sp.insert(p);
            ad.push_back(AliceData[p]), bd.push_back(BobData[p]);
        }
        for(int i = cB; i < cA; i++) ad.push_back(AliceData[i]);
        for(int i = cA; i < cB; i++) bd.push_back(BobData[i]);
        sort(ad.begin(), ad.end()), sort(bd.begin(), bd.end());
        shuffle(AliceData.begin(), AliceData.end(), rng);
        shuffle(BobData.begin(), BobData.end(), rng);
    // }

    // AliceData = {233613855, 84053574, 955552696, 443195438};
    // BobData = {955552696, 233613855, 84053574, 356042544, 443195438};
    // print(AliceData), print(BobData);

    auto check = [&](info &v) {
        return 0 <= v.sz && v.sz <= 2 * D && v.p.size() <= D;
    } ;
    double EncodeBegin = clock();
        info Alice = Encode(AliceData);
    double EncodeEnd = clock();
        info Bob = Encode(BobData);
        assert(check(Alice) && check(Bob));
    double DecodeBegin = clock();
        auto DiffRes = Decode(Alice, Bob);
    double DecodeEnd = clock();

    if(DiffRes.index() == 1) puts("ERROR"), exit(0);
    auto Diff = get<pair<vector<int>, vector<int> > >(DiffRes);
    auto ansad = Diff.first, ansbd = Diff.second;
    assert(ansad.size() == ad.size() && ansbd.size() == bd.size());
    for(int i = 0; i < (int)ad.size(); i++) assert(ad[i] == ansad[i]);
    for(int i = 0; i < (int)ad.size(); i++) assert(bd[i] == ansbd[i]);
    puts("AC");
    // puts("Alice - Bob : ");
    // for(int x : Diff.first) printf("%d ", x);
    // puts("");
    // puts("Bob - Alice : ");
    // for(int x : Diff.second) printf("%d ", x);
    // puts("");

    // printf("Running Time of Encode : %.3lf s\n", (EncodeEnd - EncodeBegin) / CLOCKS_PER_SEC);
    // printf("Running Time of Decode : %.3lf s\n", (DecodeEnd - DecodeBegin) / CLOCKS_PER_SEC);

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