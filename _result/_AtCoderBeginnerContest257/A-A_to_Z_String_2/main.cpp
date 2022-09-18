#include <bits/stdc++.h>
#define out(X) cout << (X) << endl;
#ifdef __LOCAL
#define DBG(X) cout << #X << " = " << (X) << endl;
#else
#define DBG(X)
#endif

using namespace std;

int main() {
#ifdef __LOCAL
    freopen("input", "r", stdin);
#endif

    int n, x;

    cin >> n;
    cin >> x;

    char A = 'A';

    out((x > n) ? char(A+(x-1)/n) : A);
}
