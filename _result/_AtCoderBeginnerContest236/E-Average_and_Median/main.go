// E - Average and Median
// https://atcoder.jp/contests/abc236/tasks/abc236_e
package main

import (
	"bufio"
	"container/heap"
	"fmt"
	"io/ioutil"
	"math"
	"math/bits"
	"os"
	"sort"
	"strconv"
	"strings"
)

var sc = bufio.NewScanner(os.Stdin)
var wtr = bufio.NewWriter(os.Stdout)

func main() {

	defer flush()

	n := ni()
	as := nis(n)
	minA := inf
	maxA := 0

	a3s := make([]int, n)
	for i, a := range as {
		a3s[i] = a * 1000
		minA = min(minA, a)
		maxA = max(maxA, a)
	}

	ans1 := bs(minA*1000, maxA*1000, func(c int) bool {
		return poss1(a3s, n, c)
	})
	out(float64(ans1) / 1000)
	ans2 := bs(minA, maxA, func(c int) bool {
		return poss2(as, n, c)
	})
	out(ans2)
}

func poss1(as []int, n, c int) bool {
	dp := make([][2]int, n+1)
	for i := 0; i < n+1; i++ {
		dp[i] = [2]int{-inf, -inf}
	}
	dp[1][1] = as[0] - c
	dp[1][0] = 0
	for i := 1; i < n; i++ {
		dp[i+1][1] = max(dp[i+1][1], dp[i][0]+as[i]-c)
		dp[i+1][1] = max(dp[i+1][1], dp[i][1]+as[i]-c)
		dp[i+1][0] = max(dp[i+1][0], dp[i][1])
	}
	return dp[n][0] >= 0 || dp[n][1] >= 0
}

func poss2(as []int, n, c int) bool {
	dp := make([][2]int, n+1)
	for i := 0; i < n+1; i++ {
		dp[i] = [2]int{-inf, -inf}
	}
	if as[0] >= c {
		dp[1][1] = 1
	} else {
		dp[1][1] = -1
	}
	dp[1][0] = 0
	for i := 1; i < n; i++ {
		if as[i] >= c {
			dp[i+1][1] = max(dp[i+1][1], dp[i][0]+1)
			dp[i+1][1] = max(dp[i+1][1], dp[i][1]+1)
		} else {
			dp[i+1][1] = max(dp[i+1][1], dp[i][0]-1)
			dp[i+1][1] = max(dp[i+1][1], dp[i][1]-1)
		}

		dp[i+1][0] = max(dp[i+1][0], dp[i][1])
	}
	return dp[n][0] > 0 || dp[n][1] > 0
}

// ==================================================
// init
// ==================================================

const inf = math.MaxInt64
const mod1000000007 = 1000000007
const mod998244353 = 998244353
const mod = mod1000000007

func init() {
	sc.Buffer([]byte{}, math.MaxInt64)
	sc.Split(bufio.ScanWords)
	if len(os.Args) > 1 && os.Args[1] == "i" {
		b, e := ioutil.ReadFile("./input")
		if e != nil {
			panic(e)
		}
		sc = bufio.NewScanner(strings.NewReader(strings.Replace(string(b), " ", "\n", -1)))
	}
}

// ==================================================
// io
// ==================================================

func ni() int {
	sc.Scan()
	i, e := strconv.Atoi(sc.Text())
	if e != nil {
		panic(e)
	}
	return i
}

func ni2() (int, int) {
	return ni(), ni()
}

func ni3() (int, int, int) {
	return ni(), ni(), ni()
}

func ni4() (int, int, int, int) {
	return ni(), ni(), ni(), ni()
}

func nis(arg ...int) []int {
	n := arg[0]
	t := 0
	if len(arg) == 2 {
		t = arg[1]
	}

	a := make([]int, n)
	for i := 0; i < n; i++ {
		a[i] = ni() - t
	}
	return a
}

func ni2s(n int) ([]int, []int) {
	a := make([]int, n)
	b := make([]int, n)
	for i := 0; i < n; i++ {
		a[i], b[i] = ni2()
	}
	return a, b
}

func ni3s(n int) ([]int, []int, []int) {
	a := make([]int, n)
	b := make([]int, n)
	c := make([]int, n)
	for i := 0; i < n; i++ {
		a[i], b[i], c[i] = ni3()
	}
	return a, b, c
}

func ni4s(n int) ([]int, []int, []int, []int) {
	a := make([]int, n)
	b := make([]int, n)
	c := make([]int, n)
	d := make([]int, n)
	for i := 0; i < n; i++ {
		a[i], b[i], c[i], d[i] = ni4()
	}
	return a, b, c, d
}

func ni2a(n int) [][2]int {
	a := make([][2]int, n)
	for i := 0; i < n; i++ {
		a[i][0], a[i][1] = ni2()
	}
	return a
}

func ni3a(n int) [][3]int {
	a := make([][3]int, n)
	for i := 0; i < n; i++ {
		a[i][0], a[i][1], a[i][2] = ni3()
	}
	return a
}

func ni4a(n int) [][4]int {
	a := make([][4]int, n)
	for i := 0; i < n; i++ {
		a[i][0], a[i][1], a[i][2], a[i][3] = ni4()
	}
	return a
}

func nf() float64 {
	sc.Scan()
	f, e := strconv.ParseFloat(sc.Text(), 64)
	if e != nil {
		panic(e)
	}
	return f
}

func ns() string {
	sc.Scan()
	return sc.Text()
}

func nsis() []int {
	sc.Scan()
	s := sc.Text()
	return stois(s, '0')
}

func scani() int {
	var i int
	fmt.Scanf("%i", &i)
	return i
}

func scans() string {
	var s string
	fmt.Scanf("%s", &s)
	return s
}

func out(v ...interface{}) {
	_, e := fmt.Fprintln(wtr, v...)
	if e != nil {
		panic(e)
	}
}

func outf(f string, v ...interface{}) {
	out(fmt.Sprintf(f, v...))
}

func outwoln(v ...interface{}) {
	_, e := fmt.Fprint(wtr, v...)
	if e != nil {
		panic(e)
	}
}

func outis(sl []int) {
	r := make([]string, len(sl))
	for i, v := range sl {
		r[i] = itoa(v)
	}
	out(strings.Join(r, " "))
}

func flush() {
	e := wtr.Flush()
	if e != nil {
		panic(e)
	}
}

func nftoi(decimalLen int) int {
	sc.Scan()
	s := sc.Text()

	r := 0
	minus := strings.Split(s, "-")
	isMinus := false
	if len(minus) == 2 {
		s = minus[1]
		isMinus = true
	}

	t := strings.Split(s, ".")
	i := atoi(t[0])
	r += i * pow(10, decimalLen)
	if len(t) > 1 {
		i = atoi(t[1])
		i *= pow(10, decimalLen-len(t[1]))
		r += i
	}
	if isMinus {
		return -r
	}
	return r
}

func atoi(s string) int {
	i, e := strconv.Atoi(s)
	if e != nil {
		panic(e)
	}
	return i
}

func itoa(i int) string {
	return strconv.Itoa(i)
}

func btoi(b byte) int {
	return atoi(string(b))
}

// ==================================================
// num
// ==================================================

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func maxf(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

func maxs(a *int, b int) {
	if *a < b {
		*a = b
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func mins(a *int, b int) {
	if *a > b {
		*a = b
	}
}

func abs(a int) int {
	if a > 0 {
		return a
	}
	return -a
}

func pow(a, b int) int {
	return int(math.Pow(float64(a), float64(b)))
}

func pow2(a int) int {
	return int(math.Pow(2, float64(a)))
}

func pow10(a int) int {
	return int(math.Pow(10, float64(a)))
}

func sqrt(i int) int {
	return int(math.Sqrt(float64(i)))
}

func sqrtf(i int) float64 {
	return math.Sqrt(float64(i))
}

func ch(cond bool, ok, ng int) int {
	if cond {
		return ok
	}
	return ng
}

func mul(a, b int) (int, int) {
	if a < 0 {
		a, b = -a, -b
	}
	if a == 0 || b == 0 {
		return 0, 0
	} else if a > 0 && b > 0 && a > math.MaxInt64/b {
		return 0, +1
	} else if a < math.MinInt64/b {
		return 0, -1
	}
	return a * b, 0
}

func getAngle(x, y float64) float64 {
	return math.Atan2(y, x) * 180 / math.Pi
}

func permutation(n int, k int) int {
	if k > n || k <= 0 {
		panic(fmt.Sprintf("invalid param n:%v k:%v", n, k))
	}
	v := 1
	for i := 0; i < k; i++ {
		v *= (n - i)
	}
	return v
}

/*
for {

		// Do something

		if !nextPermutation(sort.IntSlice(x)) {
			break
		}
	}
*/
func nextPermutation(x sort.Interface) bool {
	n := x.Len() - 1
	if n < 1 {
		return false
	}
	j := n - 1
	for ; !x.Less(j, j+1); j-- {
		if j == 0 {
			return false
		}
	}
	l := n
	for !x.Less(j, l) {
		l--
	}
	x.Swap(j, l)
	for k, l := j+1, n; k < l; {
		x.Swap(k, l)
		k++
		l--
	}
	return true
}

type combFactorial struct {
	fac    []int
	facinv []int
}

func newcombFactorial(n int) *combFactorial {

	fac := make([]int, n)
	facinv := make([]int, n)
	fac[0] = 1
	facinv[0] = minvfermat(1, mod)

	for i := 1; i < n; i++ {
		fac[i] = mmul(i, fac[i-1])
		facinv[i] = minvfermat(fac[i], mod)
	}

	return &combFactorial{
		fac:    fac,
		facinv: facinv,
	}
}

func (c *combFactorial) factorial(n int) int {
	return c.fac[n]
}

func (c *combFactorial) combination(n, r int) int {
	if r > n {
		return 0
	}
	return mmul(mmul(c.fac[n], c.facinv[r]), c.facinv[n-r])
}

func (c *combFactorial) permutation(n, r int) int {
	if r > n {
		return 0
	}
	return mmul(c.fac[n], c.facinv[n-r])
}

func (c *combFactorial) homogeousProduct(n, r int) int {
	return c.combination(n-1+r, r)
}

func gcd(a, b int) int {
	if b == 0 {
		return a
	}
	return gcd(b, a%b)
}

func divisor(n int) ([]int, map[int]int) {
	sqrtn := int(math.Sqrt(float64(n)))
	c := 2
	divisor := []int{}
	divisorm := make(map[int]int)
	for {
		if n%2 != 0 {
			break
		}
		divisor = append(divisor, 2)
		divisorm[2]++
		n /= 2
	}
	c = 3
	for {
		if n%c == 0 {
			divisor = append(divisor, c)
			divisorm[c]++
			n /= c
		} else {
			c += 2
			if c > sqrtn {
				break
			}
		}
	}
	if n != 1 {
		divisor = append(divisor, n)
		divisorm[n]++
	}
	return divisor, divisorm
}

type binom struct {
	fac  []int
	finv []int
	inv  []int
}

func newbinom(n int) *binom {
	b := &binom{
		fac:  make([]int, n),
		finv: make([]int, n),
		inv:  make([]int, n),
	}
	b.fac[0] = 1
	b.fac[1] = 1
	b.inv[1] = 1
	b.finv[0] = 1
	b.finv[1] = 1
	for i := 2; i < n; i++ {
		b.fac[i] = b.fac[i-1] * i % mod
		b.inv[i] = mod - mod/i*b.inv[mod%i]%mod
		b.finv[i] = b.finv[i-1] * b.inv[i] % mod
	}
	return b
}

func (b *binom) get(n, r int) int {
	if n < r || n < 0 || r < 0 {
		return 0
	}
	return b.fac[n] * b.finv[r] % mod * b.finv[n-r] % mod
}

func matPow(a [][]int, n int) [][]int {
	r := make([][]int, len(a))
	for i := 0; i < len(a); i++ {
		r[i] = is(len(a), 0)
		r[i][i] = 1
	}
	for n > 0 {
		if n&1 != 0 {
			r = matMul(a, r)
		}
		a = matMul(a, a)
		n = n >> 1
	}
	return r
}

func matMul(a, b [][]int) [][]int {
	r := make([][]int, len(a))
	for i := 0; i < len(a); i++ {
		r[i] = is(len(b[0]), 0)
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b[0]); j++ {
			for k := 0; k < len(b); k++ {
				r[i][j] = madd(r[i][j], mmul(a[i][k], b[k][j]))
			}
		}
	}
	return r
}

// ==================================================
// mod
// ==================================================

func madd(a, b int) int {
	a += b
	if a >= mod || a <= -mod {
		a %= mod
	}
	if a < 0 {
		a += mod
	}
	return a
}

func mmul(a, b int) int {
	return a * b % mod
}

func mdiv(a, b int) int {
	a %= mod
	return a * minvfermat(b, mod) % mod
}

func mpow(a, n, m int) int {
	if m == 1 {
		return 0
	}
	r := 1
	for n > 0 {
		if n&1 == 1 {
			r = r * a % m
		}
		a, n = a*a%m, n>>1
	}
	return r
}

func minv(a, m int) int {
	p, x, u := m, 1, 0
	for p != 0 {
		t := a / p
		a, p = p, a-t*p
		x, u = u, x-t*u
	}
	x %= m
	if x < 0 {
		x += m
	}
	return x
}

// m only allow prime number
func minvfermat(a, m int) int {
	return mpow(a, m-2, mod)
}

// ==================================================
// binarysearch
// ==================================================

/*
	o = bs(0, len(sl)-1, func(c int) bool {
		return true
	})
*/
func bs(ok, ng int, f func(int) bool) int {
	if !f(ok) {
		return -1
	}
	if f(ng) {
		return ng
	}
	for abs(ok-ng) > 1 {
		mid := (ok + ng) / 2

		if f(mid) {
			ok = mid
		} else {
			ng = mid
		}
	}

	return ok
}

/*
	o = bsfl(0.0, 100.0, 100, func(c float64) bool {
		return true
	})
*/
func bsfl(ok, ng float64, c int, f func(float64) bool) float64 {
	for i := 0; i < c; i++ {

		mid := (ok + ng) / 2

		if f(mid) {
			ok = mid
		} else {
			ng = mid
		}
	}

	return ok
}

func bs3fl(low, high float64, c int, f func(float64) float64) float64 {

	for i := 0; i < c; i++ {
		c1 := (low*2 + high) / 3
		c2 := (low + high*2) / 3

		if f(c1) > f(c2) {
			low = c1
		} else {
			high = c2
		}
	}
	return low
}

// ==================================================
// bit
// ==================================================

func hasbit(a int, n int) bool {
	return (a>>uint(n))&1 == 1
}

func nthbit(a int, n int) int {
	return int((a >> uint(n)) & 1)
}

func popcount(a int) int {
	return bits.OnesCount(uint(a))
}

func bitlen(a int) int {
	return bits.Len(uint(a))
}

func xor(a, b bool) bool { return a != b }

func debugbit(n int) string {
	r := ""
	for i := bitlen(n) - 1; i >= 0; i-- {
		if n&(1<<i) != 0 {
			r += "1"
		} else {
			r += "0"
		}
	}
	return r
}

// ==================================================
// string
// ==================================================

func toLowerCase(s string) string {
	return strings.ToLower(s)
}

func toUpperCase(s string) string {
	return strings.ToUpper(s)
}

func isLower(b byte) bool {
	return 'a' <= b && b <= 'z'
}

func isUpper(b byte) bool {
	return 'A' <= b && b <= 'Z'
}

// ==================================================
// sort
// ==================================================

type sortOrder int

const (
	asc sortOrder = iota
	desc
)

func sorti(sl []int) {
	sort.Sort(sort.IntSlice(sl))
}

func sortir(sl []int) {
	sort.Sort(sort.Reverse(sort.IntSlice(sl)))
}

func sorts(sl []string) {
	sort.Slice(sl, func(i, j int) bool {
		return sl[i] < sl[j]
	})
}

type Sort2ArOptions struct {
	keys   []int
	orders []sortOrder
}

type Sort2ArOption func(*Sort2ArOptions)

func opt2ar(key int, order sortOrder) Sort2ArOption {
	return func(args *Sort2ArOptions) {
		args.keys = append(args.keys, key)
		args.orders = append(args.orders, order)
	}
}

// sort2ar(sl,opt2ar(1,asc))
// sort2ar(sl,opt2ar(0,asc),opt2ar(1,asc))
func sort2ar(sl [][2]int, setters ...Sort2ArOption) {
	args := &Sort2ArOptions{}

	for _, setter := range setters {
		setter(args)
	}

	sort.Slice(sl, func(i, j int) bool {
		for idx, key := range args.keys {
			if sl[i][key] == sl[j][key] {
				continue
			}
			switch args.orders[idx] {
			case asc:
				return sl[i][key] < sl[j][key]
			case desc:
				return sl[i][key] > sl[j][key]
			}
		}
		return true
	})
}

// ==================================================
// slice
// ==================================================

func is(l int, def int) []int {
	sl := make([]int, l)
	for i := 0; i < l; i++ {
		sl[i] = def
	}
	return sl
}

func i2s(l, m int, def int) [][]int {
	sl := make([][]int, l)
	for i := 0; i < l; i++ {
		sl[i] = make([]int, m)
		for j := 0; j < m; j++ {
			sl[i][j] = def
		}
	}
	return sl
}

// out(stois("abcde", 'a'))
// out(stois("abcde", 'a'-1))
// out(stois("12345", '0'))
func stois(s string, baseRune rune) []int {
	r := make([]int, len(s))
	for i, v := range s {
		r[i] = int(v - baseRune)
	}
	return r
}

func istos(s []int, baseRune rune) string {
	r := ""
	for _, v := range s {
		r += string(v + int(baseRune))
	}
	return r
}

func reverse(sl []interface{}) {
	for i, j := 0, len(sl)-1; i < j; i, j = i+1, j-1 {
		sl[i], sl[j] = sl[j], sl[i]
	}
}

func reversei(sl []int) {
	for i, j := 0, len(sl)-1; i < j; i, j = i+1, j-1 {
		sl[i], sl[j] = sl[j], sl[i]
	}
}

func uniquei(sl []int) []int {
	hist := make(map[int]struct{})
	j := 0
	rsl := make([]int, len(sl))
	for i := 0; i < len(sl); i++ {
		if _, ok := hist[sl[i]]; ok {
			continue
		}
		rsl[j] = sl[i]
		hist[sl[i]] = struct{}{}
		j++
	}
	return rsl[:j]
}

// coordinate compression
func cocom(sl []int) ([]int, map[int]int) {
	rsl := uniquei(sl)
	sorti(rsl)
	rm := make(map[int]int)
	for i := 0; i < len(rsl); i++ {
		rm[rsl[i]] = i
	}
	return rsl, rm
}

func popBack(sl []int) (int, []int) {
	return sl[len(sl)-1], sl[:len(sl)-1]
}

func addIdx(pos, v int, sl []int) []int {
	if len(sl) == pos {
		sl = append(sl, v)
		return sl
	}
	sl = append(sl[:pos+1], sl[pos:]...)
	sl[pos] = v
	return sl
}

func delIdx(pos int, sl []int) []int {
	return append(sl[:pos], sl[pos+1:]...)
}

// find x of sl[x] < v. return -1 if no lowerbound found
func lowerBound(v int, sl []int) int {
	if len(sl) == 0 {
		panic("slise len is zero")
	}
	idx := bs(0, len(sl)-1, func(c int) bool {
		return sl[c] < v
	})
	return idx
}

// find x of v < sl[x]. return len(sl) if no upperbound found
func upperBound(v int, sl []int) int {
	if len(sl) == 0 {
		panic("slise len is zero")
	}
	idx := bs(0, len(sl)-1, func(c int) bool {
		return sl[c] <= v
	})
	return idx + 1
}

// ==================================================
// point
// ==================================================

type point struct {
	x int
	y int
}

type pointf struct {
	x float64
	y float64
}

func (p point) isValid(x, y int) bool {
	return 0 <= p.x && p.x < x && 0 <= p.y && p.y < y
}

func pointAdd(a, b point) point {
	return point{x: a.x + b.x, y: a.y + b.y}
}

func pointDist(a, b point) float64 {
	return sqrtf((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y))
}

func pointDistDouble(a, b point) int {
	return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)
}

func pointfDist(a, b pointf) float64 {
	return math.Sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y))
}

// ==================================================
// queue
// ==================================================

/*
	q := list.New()
	q.PushBack(val)
	e := q.Front()
	for e != nil {
		t := e.Value.(int)

		// Do something

		e = e.Next()
    }
*/

// ==================================================
// heap
// ==================================================

/*
ih := newIntHeap(asc)
ih.Push(v)

	for !ih.IsEmpty() {
		v := ih.Pop(h)
	}
*/
type IntHeap struct {
	sum int
	pq  *pq
}

func newIntHeap(order sortOrder) *IntHeap {
	ih := &IntHeap{}
	ih.pq = newpq([]compFunc{func(p, q interface{}) int {
		if p.(int) == q.(int) {
			return 0
		}
		if order == asc {
			if p.(int) < q.(int) {
				return -1
			} else {
				return 1
			}
		} else {
			if p.(int) > q.(int) {
				return -1
			} else {
				return 1
			}
		}
	}})
	heap.Init(ih.pq)
	return ih
}
func (ih *IntHeap) Push(x int) {
	ih.sum += x
	heap.Push(ih.pq, x)
}

func (ih *IntHeap) Pop() int {
	v := heap.Pop(ih.pq).(int)
	ih.sum -= v
	return v
}

func (ih *IntHeap) Len() int {
	return ih.pq.Len()
}

func (ih *IntHeap) IsEmpty() bool {
	return ih.pq.Len() == 0
}

func (ih *IntHeap) GetRoot() int {
	return ih.pq.GetRoot().(int)
}

func (ih *IntHeap) GetSum() int {
	return ih.sum
}

/*
h := &OrgIntHeap{}
heap.Init(h)

heap.Push(h, v)

	for !h.IsEmpty() {
		v = heap.Pop(h).(int)
	}
*/
type OrgIntHeap []int

func (h OrgIntHeap) Len() int { return len(h) }

// get from bigger
// func (h OrgIntHeap) Less(i, j int) bool { return h[i] > h[j] }
func (h OrgIntHeap) Less(i, j int) bool { return h[i] < h[j] }
func (h OrgIntHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }

func (h *OrgIntHeap) Push(x interface{}) {
	*h = append(*h, x.(int))
}

func (h *OrgIntHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func (h *OrgIntHeap) IsEmpty() bool {
	return h.Len() == 0
}

// h.Min().(int)
func (h *OrgIntHeap) Min() interface{} {
	return (*h)[0]
}

/*
	type pqst struct {
		x int
		y int
	}

	pq := newpq([]compFunc{func(p, q interface{}) int {
		if p.(pqst).x != q.(pqst).x {
			// get from bigger
			// if p.(pqst).x > q.(pqst).x {
			if p.(pqst).x < q.(pqst).x {
				return -1
			} else {
				return 1
			}
		}
		if p.(pqst).y != q.(pqst).y {
			// get from bigger
			// if p.(pqst).y > q.(pqst).y {
			if p.(pqst).y < q.(pqst).y {
				return -1
			} else {
				return 1
			}
		}
		return 0
	}})
	heap.Init(pq)
	heap.Push(pq, pqst{x: 1, y: 1})
	for !pq.IsEmpty() {
		v := heap.Pop(pq).(pqst)
	}
*/

type pq struct {
	arr   []interface{}
	comps []compFunc
}

type compFunc func(p, q interface{}) int

func newpq(comps []compFunc) *pq {
	return &pq{
		comps: comps,
	}
}

func (pq pq) Len() int {
	return len(pq.arr)
}

func (pq pq) Swap(i, j int) {
	pq.arr[i], pq.arr[j] = pq.arr[j], pq.arr[i]
}

func (pq pq) Less(i, j int) bool {
	for _, comp := range pq.comps {
		result := comp(pq.arr[i], pq.arr[j])
		switch result {
		case -1:
			return true
		case 1:
			return false
		case 0:
			continue
		}
	}
	return true
}

func (pq *pq) Push(x interface{}) {
	pq.arr = append(pq.arr, x)
}

func (pq *pq) Pop() interface{} {
	n := pq.Len()
	item := pq.arr[n-1]
	pq.arr = pq.arr[:n-1]
	return item
}

func (pq *pq) IsEmpty() bool {
	return pq.Len() == 0
}

// pq.GetRoot().(edge)
func (pq *pq) GetRoot() interface{} {
	return pq.arr[0]
}

// ==================================================
// cusum
// ==================================================

type cusum struct {
	l int
	s []int
}

func newcusum(sl []int) *cusum {
	c := &cusum{}
	c.l = len(sl)
	c.s = make([]int, len(sl)+1)
	for i, v := range sl {
		c.s[i+1] = c.s[i] + v
	}
	return c
}

// get sum f <= i && i <= t
func (c *cusum) getRange(f, t int) int {
	if f > t || f >= c.l {
		return 0
	}
	return c.s[t+1] - c.s[f]
}

// get sum 0 to i
func (c *cusum) get(i int) int {
	return c.s[i+1]
}

func (c *cusum) upperBound(i int) int {
	return upperBound(i, c.s)
}

func (c *cusum) lowerBound(i int) int {
	return lowerBound(i, c.s)
}

/*
	mp := make([][]int, n)
	for i := 0; i < k; i++ {
		mp[i] = make([]int, m)
	}
	cusum2d := newcusum2d(sl)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			t:=cusum2d.get(0, 0, i, j)
		}
	}
*/

type cusum2d struct {
	s [][]int
}

func newcusum2d(sl [][]int) *cusum2d {
	c := &cusum2d{}
	n := len(sl)
	m := len(sl[0])
	c.s = make([][]int, n+1)
	for i := 0; i < n+1; i++ {
		c.s[i] = make([]int, m+1)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			c.s[i+1][j+1] = c.s[i+1][j] + c.s[i][j+1] - c.s[i][j]
			c.s[i+1][j+1] += sl[i][j]
		}
	}
	return c
}

// x1 <= x <= x2, y1 <= y <= y2
func (c *cusum2d) get(x1, y1, x2, y2 int) int {
	return c.s[x2][y2] + c.s[x1][y1] - c.s[x1][y2] - c.s[x2][y1]
}

// ==================================================
// union find
// ==================================================

type unionFind struct {
	par []int
}

func newUnionFind(n int) *unionFind {
	u := &unionFind{
		par: make([]int, n),
	}
	for i := range u.par {
		u.par[i] = -1
	}
	return u
}

func (u *unionFind) root(x int) int {
	if u.par[x] < 0 {
		return x
	}
	u.par[x] = u.root(u.par[x])
	return u.par[x]
}

func (u *unionFind) unite(x, y int) {
	x = u.root(x)
	y = u.root(y)
	if x == y {
		return
	}
	if u.size(x) < u.size(y) {
		x, y = y, x
	}
	u.par[x] += u.par[y]
	u.par[y] = x
}

func (u *unionFind) issame(x, y int) bool {
	if u.root(x) == u.root(y) {
		return true
	}
	return false
}

func (u *unionFind) size(x int) int {
	return -u.par[u.root(x)]
}

// ==================================================
// bit
// ==================================================

type bit struct {
	n int
	b []int
}

func newbit(n int) *bit {
	return &bit{
		n: n + 1,
		b: make([]int, n+1),
	}
}

func (b *bit) culc(i, j int) int {
	return i + j
	//return madd(i, j)
}

func (b *bit) add(i, x int) {
	for i++; i < b.n && i > 0; i += i & -i {
		b.b[i] = b.culc(b.b[i], x)
	}
}

func (b *bit) sum(i int) int {
	ret := 0
	for i++; i > 0; i -= i & -i {
		ret = b.culc(ret, b.b[i])
	}
	return ret
}

// l <= x < r
func (b *bit) rangesum(l, r int) int {
	return b.culc(b.sum(r-1), -b.sum(l-1))
}

func (b *bit) lowerBound(x int) int {
	idx, k := 0, 1
	for k < b.n {
		k <<= 1
	}
	for k >>= 1; k > 0; k >>= 1 {
		if idx+k < b.n && b.b[idx+k] < x {
			x -= b.b[idx+k]
			idx += k
		}
	}
	return idx
}

// K:キーの型, V:値の型
func NewRBMAP() *RBMAP {
	return &RBMAP{}
}

type RBMAPColor int

///////////////////////////////////////////////////////////////////////////
// 共通定義
///////////////////////////////////////////////////////////////////////////

// R:赤, B:黒, Error:debug 用
const (
	RBMAPColorR RBMAPColor = iota
	RBMAPColorB
	RBMAPColorError
)

type RBMAPNode struct { // ノードの型
	color RBMAPColor // そのノードの色
	key   int        // そのノードのキー
	value int        // そのノードの値
	lst   *RBMAPNode // 左部分木
	rst   *RBMAPNode // 右部分木
}

func NewNode(color RBMAPColor, key int, value int) *RBMAPNode {
	return &RBMAPNode{
		color: color,
		key:   key,
		value: value,
	}
}

type RBMAP struct {
	root   *RBMAPNode // 赤黒木の根
	change bool       // 修正が必要かを示すフラグ(true:必要, false:不要)
	lmax   int        // 左部分木のキーの最大値
	value  int        // lmax に対応する値
}

// ノード n が赤かチェックする
func (n *RBMAPNode) isR() bool {
	return n != nil && n.color == RBMAPColorR
}

// ノード n が黒かチェックする
func (n *RBMAPNode) isB() bool {
	return n != nil && n.color == RBMAPColorB
}

// ２分探索木 v の左回転。回転した木を返す
func rotateL(v *RBMAPNode) *RBMAPNode {
	u := v.rst
	t2 := u.lst
	u.lst = v
	v.rst = t2
	return u
}

// ２分探索木 u の右回転。回転した木を返す
func rotateR(u *RBMAPNode) *RBMAPNode {
	v := u.lst
	t2 := v.rst
	v.rst = u
	u.lst = t2
	return v
}

// ２分探索木 t の二重回転(左回転 -> 右回転)。回転した木を返す
func rotateLR(t *RBMAPNode) *RBMAPNode {
	t.lst = rotateL(t.lst)
	return rotateR(t)
}

// ２分探索木 t の二重回転(右回転 -> 左回転)。回転した木を返す
func rotateRL(t *RBMAPNode) *RBMAPNode {
	t.rst = rotateR(t.rst)
	return rotateL(t)
}

///////////////////////////////////////////////////////////////////////////
// insert(挿入)
///////////////////////////////////////////////////////////////////////////

// エントリー(key, x のペア)を挿入する
func (m *RBMAP) Insert(key int, x int) {
	m.root = m.insertSub(m.root, key, x)
	m.root.color = RBMAPColorB
}

func (m *RBMAP) Increment(key int) {
	m.Insert(key, m.Lookup(key)+1)
}

func (m *RBMAP) insertSub(t *RBMAPNode, key int, x int) *RBMAPNode {
	if t == nil {
		m.change = true
		return NewNode(RBMAPColorR, key, x)
	}
	cmp := 0
	if key > t.key {
		cmp = 1
	}
	if key < t.key {
		cmp = -1
	}
	// cmp > 0 を cmp >= 0にするとmultisetになる
	if cmp < 0 {
		t.lst = m.insertSub(t.lst, key, x)
		return m.balance(t)
	} else if cmp >= 0 {
		t.rst = m.insertSub(t.rst, key, x)
		return m.balance(t)
	}
	m.change = false
	t.value = x
	return t
}

// エントリー挿入に伴う赤黒木の修正(パターンマッチ)
func (m *RBMAP) balance(t *RBMAPNode) *RBMAPNode {
	if !m.change {
		return t
	} else if !t.isB() {
		return t // 根が黒でないなら何もしない
	} else if t.lst.isR() && t.lst.lst.isR() {
		t = rotateR(t)
		t.lst.color = RBMAPColorB
	} else if t.lst.isR() && t.lst.rst.isR() {
		t = rotateLR(t)
		t.lst.color = RBMAPColorB
	} else if t.rst.isR() && t.rst.lst.isR() {
		t = rotateRL(t)
		t.rst.color = RBMAPColorB
	} else if t.rst.isR() && t.rst.rst.isR() {
		t = rotateL(t)
		t.rst.color = RBMAPColorB
	} else {
		m.change = false
	}
	return t
}

///////////////////////////////////////////////////////////////////////////
// delete(削除)
///////////////////////////////////////////////////////////////////////////

// key で指すエントリー(ノード)を削除する
func (m *RBMAP) Delete(key int) {
	m.root = m.deleteSub(m.root, key)
	if m.root != nil {
		m.root.color = RBMAPColorB
	}
}

func (m *RBMAP) Decrement(key int) {
	count := m.Lookup(key)
	count--
	if count <= 0 {
		m.Delete(key)
	} else {
		m.Insert(key, count)
	}
}

func (m *RBMAP) deleteSub(t *RBMAPNode, key int) *RBMAPNode {
	if t == nil {
		m.change = false
		return nil
	}
	cmp := 0
	if key > t.key {
		cmp = 1
	}
	if key < t.key {
		cmp = -1
	}
	if cmp < 0 {
		t.lst = m.deleteSub(t.lst, key)
		return m.balanceL(t)
	} else if cmp > 0 {
		t.rst = m.deleteSub(t.rst, key)
		return m.balanceR(t)
	} else {
		if t.lst == nil {
			switch t.color {
			case RBMAPColorR:
				m.change = false
				break
			case RBMAPColorB:
				m.change = true
				break
			}
			return t.rst // 右部分木を昇格させる
		} else {
			t.lst = m.deleteMax(t.lst) // 左部分木の最大値を削除する
			t.key = m.lmax             // 左部分木の削除した最大値で置き換える
			t.value = m.value
			return m.balanceL(t)
		}
	}
}

// 部分木 t の最大値のノードを削除する
// 戻り値は削除により修正された部分木
// 削除した最大値を lmax に保存する
func (m *RBMAP) deleteMax(t *RBMAPNode) *RBMAPNode {
	if t.rst != nil {
		t.rst = m.deleteMax(t.rst)
		return m.balanceR(t)
	} else {
		m.lmax = t.key // 部分木のキーの最大値を保存
		m.value = t.value
		switch t.color {
		case RBMAPColorR:
			m.change = false
			break
		case RBMAPColorB:
			m.change = true
			break
		}
		return t.lst // 左部分木を昇格させる
	}
}

// 左部分木のノード削除に伴う赤黒木の修正(パターンマッチ)
// 戻り値は修正された木
func (m *RBMAP) balanceL(t *RBMAPNode) *RBMAPNode {
	if !m.change {
		return t // 修正なしと赤ノード削除の場合はここ

	} else if t.rst.isB() && t.rst.lst.isR() {
		rb := t.color
		t = rotateRL(t)
		t.color = rb
		t.lst.color = RBMAPColorB
		m.change = false
	} else if t.rst.isB() && t.rst.rst.isR() {
		rb := t.color
		t = rotateL(t)
		t.color = rb
		t.lst.color = RBMAPColorB
		t.rst.color = RBMAPColorB
		m.change = false
	} else if t.rst.isB() {
		rb := t.color
		t.color = RBMAPColorB
		t.rst.color = RBMAPColorR
		m.change = (rb == RBMAPColorB)
	} else if t.rst.isR() {
		t = rotateL(t)
		t.color = RBMAPColorB
		t.lst.color = RBMAPColorR
		t.lst = m.balanceL(t.lst)
		m.change = false
	} else { // 黒ノード削除の場合、ここはありえない
		panic("(L) This program is buggy")
	}
	return t
}

// 右部分木のノード削除に伴う赤黒木の修正(パターンマッチ)
// 戻り値は修正された木
func (m *RBMAP) balanceR(t *RBMAPNode) *RBMAPNode {
	if !m.change {
		return t // 修正なしと赤ノード削除の場合はここ
	} else if t.lst.isB() && t.lst.rst.isR() {
		rb := t.color
		t = rotateLR(t)
		t.color = rb
		t.rst.color = RBMAPColorB
		m.change = false
	} else if t.lst.isB() && t.lst.lst.isR() {
		rb := t.color
		t = rotateR(t)
		t.color = rb
		t.lst.color = RBMAPColorB
		t.rst.color = RBMAPColorB
		m.change = false
	} else if t.lst.isB() {
		rb := t.color
		t.color = RBMAPColorB
		t.lst.color = RBMAPColorR
		m.change = (rb == RBMAPColorB)
	} else if t.lst.isR() {
		t = rotateR(t)
		t.color = RBMAPColorB
		t.rst.color = RBMAPColorR
		t.rst = m.balanceR(t.rst)
		m.change = false
	} else { // 黒ノード削除の場合、ここはありえない
		panic("(R) This program is buggy")
	}
	return t
}

///////////////////////////////////////////////////////////////////////////
// member(検索)等
///////////////////////////////////////////////////////////////////////////

// キーの検索。ヒットすれば true、しなければ false
func (m *RBMAP) Member(key int) bool {
	t := m.root
	for t != nil {
		cmp := 0
		if key > t.key {
			cmp = 1
		}
		if key < t.key {
			cmp = -1
		}
		if cmp < 0 {
			t = t.lst
		} else if cmp > 0 {
			t = t.rst
		} else {
			return true
		}
	}
	return false
}

// 指定されたキーより大きいキーの検索
func (m *RBMAP) UpperBound(key int) (int, bool) {
	t := m.root
	r := 0
	hasKey := false
	for t != nil {
		cmp := 0
		if key > t.key {
			cmp = 1
		}
		if key < t.key {
			cmp = -1
		}
		if cmp < 0 {
			if !hasKey || t.key < r {
				hasKey = true
				r = t.key
			}
			t = t.lst
		} else {
			t = t.rst
		}
	}
	return r, hasKey
}

// 指定されたキーより小さいキーの検索
func (m *RBMAP) LowerBound(key int) (int, bool) {
	t := m.root
	l := 0
	hasKey := false
	for t != nil {
		cmp := 0
		if key > t.key {
			cmp = 1
		}
		if key < t.key {
			cmp = -1
		}
		if cmp > 0 {
			if !hasKey || t.key > l {
				hasKey = true
				l = t.key
			}
			t = t.rst
		} else {
			t = t.lst
		}
	}
	return l, hasKey
}

// キーから値を得る。キーがヒットしない場合は nil を返す
func (m *RBMAP) Lookup(key int) int {
	t := m.root
	for t != nil {
		cmp := 0
		if key > t.key {
			cmp = 1
		}
		if key < t.key {
			cmp = -1
		}
		if cmp < 0 {
			t = t.lst
		} else if cmp > 0 {
			t = t.rst
		} else {
			return t.value
		}
	}
	return 0
}

// マップが空なら true、空でないなら false
func (m *RBMAP) IsEmpty() bool {
	return m.root == nil
}

// マップを空にする
func (m *RBMAP) Clear() {
	m.root = nil
}

// キーのリスト
func (m *RBMAP) Keys() []int {
	al := []int{}
	al = m.keysSub(m.root, al)
	return al
}

// 値のリスト
func (m *RBMAP) Values() []int {
	al := []int{}
	al = m.valuesSub(m.root, al)
	return al
}

// マップのサイズ
func (m *RBMAP) Size() int {
	return len(m.Keys())
}

// キーの最小値
func (m *RBMAP) Min() int {
	t := m.root
	if t == nil {
		return 0
	}
	for t.lst != nil {
		t = t.lst
	}
	return t.key
}

// キーの最大値
func (m *RBMAP) Max() int {
	t := m.root
	if t == nil {
		return 0
	}
	for t.rst != nil {
		t = t.rst
	}
	return t.key
}

func (m *RBMAP) keysSub(t *RBMAPNode, al []int) []int {
	if t != nil {
		al = m.keysSub(t.lst, al)
		al = append(al, t.key)
		al = m.keysSub(t.rst, al)
	}
	return al
}

func (m *RBMAP) valuesSub(t *RBMAPNode, al []int) []int {
	if t != nil {
		al = m.valuesSub(t.lst, al)
		al = append(al, t.value)
		al = m.valuesSub(t.rst, al)
	}
	return al
}

///////////////////////////////////////////////////////////////////////////
// debug 用ルーチン
///////////////////////////////////////////////////////////////////////////

// 赤黒木をグラフ文字列に変換する
func (m *RBMAP) String() string {
	return m.toGraph("", "", m.root)
}

// 赤黒木のバランスが取れているか確認する
func (m *RBMAP) Balanced() bool {
	return m.blackHeight(m.root) >= 0
}

// 赤黒木の配色が正しいか確認する
func (m *RBMAP) ColorValid() bool {
	return m.colorChain(m.root) == RBMAPColorB
}

// ２分探索木になっているか確認する
func (m *RBMAP) BstValid() bool {
	return m.bstValidSub(m.root)
}

func (m *RBMAP) toGraph(head string, bar string, t *RBMAPNode) string {
	graph := ""
	if t != nil {
		graph += m.toGraph(head+"　　", "／", t.rst)
		node := ""
		switch t.color {
		case RBMAPColorR:
			node = "R"
		case RBMAPColorB:
			node = "B"
		}
		node += fmt.Sprintf(":%v:%v", t.key, t.value)
		graph += fmt.Sprintf("%v%v%v\n", head, bar, node)
		graph += m.toGraph(head+"　　", "＼", t.lst)
	}
	return graph
}

func (m *RBMAP) blackHeight(t *RBMAPNode) int {
	if t == nil {
		return 0
	}
	nl := m.blackHeight(t.lst)
	nr := m.blackHeight(t.rst)
	if nl < 0 || nr < 0 || nl != nr {
		return -1
	}
	if t.color == RBMAPColorB {
		return nl + 1
	}
	return nl
}

func (m *RBMAP) colorChain(t *RBMAPNode) RBMAPColor {
	if t == nil {
		return RBMAPColorB
	}
	p := t.color
	cl := m.colorChain(t.lst)
	cr := m.colorChain(t.rst)
	if cl == RBMAPColorError || cr == RBMAPColorError {
		return RBMAPColorError
	}
	if p == RBMAPColorB {
		return p
	}
	if p == RBMAPColorR && cl == RBMAPColorB && cr == RBMAPColorB {
		return p
	}
	return RBMAPColorError
}

func (m *RBMAP) bstValidSub(t *RBMAPNode) bool {
	if t == nil {
		return true
	}
	flag := m.small(t.key, t.lst) && m.large(t.key, t.rst)
	return flag && m.bstValidSub(t.lst) && m.bstValidSub(t.rst)
}

func (m *RBMAP) small(key int, t *RBMAPNode) bool {
	if t == nil {
		return true
	}
	flag := key > t.key
	return flag && m.small(key, t.lst) && m.small(key, t.rst)
}

func (m *RBMAP) large(key int, t *RBMAPNode) bool {
	if t == nil {
		return true
	}
	flag := key < t.key
	return flag && m.large(key, t.lst) && m.large(key, t.rst)
}
