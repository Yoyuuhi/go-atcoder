package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"strconv"
	"strings"
)

var sc = bufio.NewScanner(os.Stdin)
var wtr = bufio.NewWriter(os.Stdout)

type point struct {
	x int
	y int
}

func main() {

	defer flush()

	n, m := ni2()
	paths := make([]point, 3)
	for i := 0; i < 3; i++ {
		paths[i] = point{x: ni(), y: ni()}
	}

	ngs := map[int]map[int]bool{}
	for i := 0; i < m; i++ {
		x, y := ni2()
		if _, e := ngs[x]; !e {
			ngs[x] = map[int]bool{}
		}

		ngs[x][y] = true
	}

	var dp [301][301][301]int
	dp[0][0][0] = 1

	ans := 0
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			for k := 0; k <= i-j; k++ {
				tX := j*paths[0].x + k*paths[1].x + (i-j-k)*paths[2].x
				tY := j*paths[0].y + k*paths[1].y + (i-j-k)*paths[2].y
				for iP, path := range paths {
					toX := tX + path.x
					toY := tY + path.y
					if ngs[toX][toY] {
						continue
					}

					switch iP {
					case 0:
						dp[i+1][j+1][k] += dp[i][j][k]
						dp[i+1][j+1][k] %= mod998244353
					case 1:
						dp[i+1][j][k+1] += dp[i][j][k]
						dp[i+1][j][k+1] %= mod998244353
					case 2:
						dp[i+1][j][k] += dp[i][j][k]
						dp[i+1][j][k] %= mod998244353
					}

					if i+1 == n {
						ans += dp[i][j][k]
						ans %= mod998244353
					}
				}
			}
		}
	}

	out(ans)
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

func flush() {
	e := wtr.Flush()
	if e != nil {
		panic(e)
	}
}

func out(v ...interface{}) {
	_, e := fmt.Fprintln(wtr, v...)
	if e != nil {
		panic(e)
	}
}
