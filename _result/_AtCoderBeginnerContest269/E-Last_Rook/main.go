// https://atcoder.jp/contests/abc269/tasks/abc269_e
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

func main() {

	defer flush()

	n := ni()

	ansTop := simulate(2, n)
	ansLeft := simulate(1, n)

	out("!", ansTop, ansLeft)
}

func simulate(mode int, n int) int {
	left := 1
	right := n

	for left < right {
		mid := (left + right) / 2

		switch mode {
		case 1:
			fmt.Printf("? %v %v %v %v\n", 1, n, left, mid)
		case 2:
			fmt.Printf("? %v %v %v %v\n", left, mid, 1, n)
		}

		leftNum := ni()

		rightMid := mid
		if (left+right)%2 != 0 {
			rightMid++
		}

		if leftNum == mid-left+1 {
			left = rightMid
		} else {
			right = mid
		}
	}
	return left
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

func out(v ...interface{}) {
	_, e := fmt.Fprintln(wtr, v...)
	if e != nil {
		panic(e)
	}
}

func flush() {
	e := wtr.Flush()
	if e != nil {
		panic(e)
	}
}
