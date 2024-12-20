
<p align="center">
  <img src="https://github.com/user-attachments/assets/0a4dd5ab-85e7-4eaa-ad34-dd9b9b078b66" width=30%>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![C++](https://img.shields.io/badge/C++-00599C?style=flat-square&logo=C%2B%2B&logoColor=white)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/fockl/FMMA/actions.yml?branch=main)
[![test](https://github.com/fockl/FMMA/actions/workflows/actions.yml/badge.svg)](https://github.com/fockl/FMMA/actions/workflows/actions.yml)
![GitHub License](https://img.shields.io/github/license/fockl/FMMA)
[![Open in Visual Studio Code](https://img.shields.io/static/v1?logo=visualstudiocode&label=&message=Open%20in%20Visual%20Studio%20Code&labelColor=2c2c32&color=007acc&logoColor=007acc)](https://vscode.dev/github/fockl/FMMA)

[English](#fmmaenglish)

# FMMA

任意次元の変数 $x_i$, $y_j$ と任意の関数 $f$ について、

``` math

c_i = \sum_{j} w_j f(x_i, y_j)

```

を高速に計算するためのライブラリ

# インストール

## C++
cmakeを用いた場合、以下のようにしてインストール出来る
```sh
cmake -B build
cmake --build build
cmake --install build
```

BLASを用いて高速化する場合は、build時に
```sh
cmake -B build -DFMMA_USE_BLAS=ON
```
とする

## Python
pip を用いてインストール可能

```sh
pip install pyfmma
```

もしくは

```sh
pip install git+https://github.com/fockl/FMMA.git
```

cmakeを用いてより詳しく条件を設定したい場合は

```sh
cmake -B build
cmake --build build
```

をした後 python 側で

```python
import build.pyfmma
```

をする

# 使い方

C++の場合、
```c++
fmma::FMMA<double, 3> fmma;
fmma.set_fn([](auto x, auto y){return 1.0/(x[0]-y[0]);}); 
fmma.set_solver_type("fmm");
fmma.solve(target, source_weight, source, ans);
```

のようにして使用する

詳しくはtutorial参照

現在はsolverとして`exact`, `nrnmm`, `tree`, `fmm`が実装済み

$O(n(x)) = O(n(y)) = O(N)$の時の計算量は以下の通り：

|type|computatoin cost|
|---|---|
|exact|$O(N^2)$|
|nrnmm|$O(N\sqrt{N})$|
|tree|$O(N\log{N})$|
|fmm|$O(N)$|

# ベンチマーク結果

github-actions を用いたベンチマーク結果：

1次元の場合:

![ time ](benchmark/results/time_1.png)

![ error ](benchmark/results/error_1.png)

2次元の場合:

![ time ](benchmark/results/time_2.png)

![ error ](benchmark/results/error_2.png)

# 参考文献

- W. Fong and E. Darve. The black-box fast multipole method. Journal of Computational Physics, 228 (2009).

# FMMA(English)

FMMA is a library to calculate fastly

``` math

c_i = \sum_{j} w_j f(x_i, y_j)

```

for arbitrary function $f$ and variables $x_i$, $y_j$ in arbitrary dimension.

Benchmark results using github-actions are follows :

1D:

![ time ](benchmark/results/time_1.png)

![ error ](benchmark/results/error_1.png)

2D:

![ time ](benchmark/results/time_2.png)

![ error ](benchmark/results/error_2.png)

# Install(English)

## C++

You can install this library as follows if cmake is used:
```sh
cmake -B build
cmake --build build
cmake --install build
```

If BLAS is required,  define an argument like:
```sh
cmake -B build -DFMMA_USE_BLAS=ON
```

## Python

You can install via pip

```sh
pip install pyfmma
```

or

```sh
pip install git+https://github.com/fockl/FMMA.git
```

If you want to set details with using cmake, 

```sh
cmake -B build
cmake --build build
```

and in python

```python
import build.pyfmma
```

# Usage(English)

In C++, you can use FMMA as

```c++
fmma::FMMA<double, 3> fmma;
fmma.set_fn([](auto x, auto y){return 1.0/(x[0]-y[0]);}); 
fmma.set_solver_type("fmm");
fmma.solve(target, source_weight, source, ans);
```

For more details, see tutorials

`exact`, `nrnmm`, `tree`, `fmm` are now implemented as solver.

when $O(n(x)) = O(n(y)) = O(N)$, the computational cost are as follows:

|type|computatoin cost|
|---|---|
|exact|$O(N^2)$|
|nrnmm|$O(N\sqrt{N})$|
|tree|$O(N\log{N})$|
|fmm|$O(N)$|

# Benchmark results

Benchmark results using github-actions are as follows:

1D:

![ time ](benchmark/results/time_1.png)

![ error ](benchmark/results/error_1.png)

2D:

![ time ](benchmark/results/time_2.png)

![ error ](benchmark/results/error_2.png)

# References

- W. Fong and E. Darve. The black-box fast multipole method. Journal of Computational Physics, 228 (2009).
