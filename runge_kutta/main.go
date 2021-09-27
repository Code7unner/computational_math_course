package main

import (
	"math"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

var (
	// Cauchy problem initial conditions
	x0 = 0
	y0 = 1.0

	// X approximation array
	t = linspace(0, 4, 3)

	// Step
	h = math.Abs(t[0] - t[1])
)

func linspace(start, stop, step float64) []float64 {
	N := int(math.Ceil((stop - start) / step))
	ranger := make([]float64, N)
	for x := range ranger {
		ranger[x] = start + step*float64(x)
	}
	return ranger
}

// y' = x + y
func f(x, y float64) float64 {
	return x + y
}

// y = 2 * e^x - x - 1
func y(x float64) float64 {
	return 2.0 * math.Exp(x) - x - 1.0
}

func initSolArray() []float64 {
	result := make([]float64, len(t))
	result[x0] = y0

	return result
}

func findTrueSol() []float64 {
	result := initSolArray()

	for i := 0; i < len(t); i++ {
		result[i] = y(t[i])
	}

	return result
}

func findRungeKutta1() []float64 {
	result := initSolArray()

	for i := 0; i < len(t) - 1; i++ {
		result[i + 1] = result[i] + f(t[i], result[i]) * (t[i + 1] - t[i])
	}

	return result
}

func findRungeKutta2() []float64 {
	result := initSolArray()

	for i := 0; i < len(t) - 1; i++ {
		k1 := h * f(t[i], result[i]) * 0.5

		result[i + 1] = result[i] + h * f(result[i] + k1, t[i] + h * 0.5)
	}

	return result
}

func findRungeKutta3() []float64 {
	result := initSolArray()

	for i := 0; i < len(t) - 1; i++ {
		k1 := h * f(t[i], result[i])
		k2 := h * f(t[i] + h / 2, result[i] + k1 / 2)
		k3 := h * f(t[i] + h, result[i] - k1 + 2 * k2)

		result[i + 1] = result[i] + 1.0 / 6.0 * (k1 + 4.0 * k2 + k3)
	}

	return result
}

func findRungeKutta4() []float64 {
	result := initSolArray()

	for i := 0; i < len(t) - 1; i++ {
		k1 := h * (f(t[i], result[i]))
		k2 := h * (f((t[i] + h / 2), (result[i] + k1 / 2)))
		k3 := h * (f((t[i] + h / 2), (result[i] + k2 / 2)))
		k4 := h * (f((t[i] + h), (result[i] + k3)))
		k := (k1 + 2 * k2 + 2 * k3 + k4) / 6

		result[i + 1] = result[i] + k
	}

	return result
}

func convertToPoints(xArr, yArr []float64) (pts plotter.XYs) {
	pts = make(plotter.XYs, len(xArr))

	for i, x := range xArr {
		pts[i].X = x
		pts[i].Y = yArr[i]
	}
	return pts
}

func printPlot() error {
	p := plot.New()

	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	p.Title.Text = "Solution of y'= x + y, y(0)=1"

	if err := plotutil.AddLinePoints(p,
		"True", convertToPoints(t, findTrueSol()),
		"Runge-Kutta-1th", convertToPoints(t, findRungeKutta1()),
		"Runge-Kutta-2th", convertToPoints(t, findRungeKutta2()),
		"Runge-Kutta-3th", convertToPoints(t, findRungeKutta3()),
		"Runge-Kutta-4th", convertToPoints(t, findRungeKutta4()));
	err != nil {
		return err
	}

	if err := p.Save(4*vg.Inch, 4*vg.Inch, "runge-kutta.png"); err != nil {
		return err
	}

	return nil
}

func main()  {
	if err := printPlot(); err != nil {
		panic(err)
	}
}
