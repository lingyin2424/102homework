# GAMES 102 作业一——曲线部分

众所周知，GAMES 102 的“无境”作业框架非常难配，很多小伙伴（比如我）折腾了很久也没配好，望而却步，遂放弃 102，非常可惜。

有一日，我痛定思痛，决定要自己写一个作业框架，并分享给大家，我边学边写，遂有了这个小项目。

为了节省大家时间，像 GAMES 101 的作业一样，我只把最重要的一些功能（函数）给空了出来，大家填空就行。

该项目的内容和 102 课程的作业有些区别，仅供大家参考。

### 环境

C++ 17，Qt 5，Eigen。

windows 的小伙伴建议使用 WSL 来配置环境，不会的非常建议直接和 chatgpt 对线。

### 交互说明：

1. 在屏幕上左键点击空白处会增加一个顶点。

2. 点击一个已有的点会删除这个点。

3. 顶点可以拖拽。

### code 1

##### work 1（估计难度 30）

黄色线条

平面上有一些点，用拉格朗日插值去插这些点：

```c++
// 补全 Poly 类中的 lagrangeInterpolation 函数
struct Poly{
  	vector<double > x, y;
	double lagrangeInterpolation(double _x) const {
        assert(x.size() == y.size());
        ////////////////////////////////////////////////////// 
        // 用 this -> x 和 this -> y 对 _x 做拉格朗日插值
        
    }
}
```

##### work 2（估计难度 80）

绿色线条

使用 Gauss 基函数的线性组合 $f(x)=\sum_{i=0}^{n-1}b_i g_i(x)$  插值 $\{P_j\}$，其中：
$$
g_i(x)=\exp\left(-\frac{(x-x_i)^2}{2\sigma^2}\right)
$$
即对称轴在插值点上，$i=1,\dots,n$，默认设 $\sigma =1$ 。

**特殊交互说明：**按键盘的 “up” 键让方差增加 ，按键盘的 “down” 键让方差乘减小。

```c++
// 补全 Poly 类中的 guassInterpolation 函数
struct Poly{
  	vector<double > x, y;
    vector<double > guassInterpolation(double sigma = 1){
        ////////////////////////////////////////////////////// 
        // n 个方差为 sigma, ui = xi 的高斯函数
        // 用 this -> x 和 this -> y 对这些高斯函数插值，返回每个高斯函数的系数，即公式中的 b

    }
}
```

### code 2

##### work 1（估计难度 70）

固定幂基函数的最高次数 m (m<n)，使用最小二乘法：$\min E$，其中 $E(x)=\sum_{i=0}^{n}(y_i-f(x_i))^2$ 拟合 $\{P_j\}$。

岭回归（Ridge Regression）：对最小二乘法误差函数增加 $E_1$ 正则项，参数 $\lambda$，$\min (E+\lambda E_1)$，其中 $E_1=\sum_{i=1}^n\alpha_i^2$ 

```c++
// 补全 Poly 类中的 polynomialFit 函数
struct Poly{  
    vector<double > x, y;
    vector<double> polynomialFit(int m, double lamda = 0){
    	////////////////////////////////////////////////////// 
    	// 用 this -> x 和 this -> y 拟合一个 m 次多项式
    	// lamda 是岭回归系数，你可以先假设 lamda = 0，即没有岭回归
        // 返回拟合的多项式的系数，小的下标对应低次的系数
	}
}
```

**特殊交互说明：**按键盘的 “up” 键让 $\lambda$ 增加 ，按键盘的 “down” 键让 $\lambda$ 乘减小。按 ”left“键让拟合的点减少一个，按 ”right“键让拟合的点增加一个。

绿色线条是无岭回归的线，黄色是有岭回归的线。

### code 3

##### work 1（估计难度 0）

对曲线做参数化：

```c++
// 补全 Poly 类中的 polynomialFit 函数
struct Parameterization{
  	static vector<float > distance(const vecf& x, const vecf& y){
    	////////////////////////////////////////////////////// 
    	// 按照欧氏距离参数化 [0, 1]
	
    	return t;
  	}
  	static vector<float > average(const vecf& x, const vecf& y){
        ////////////////////////////////////////////////////// 
        // 平均参数化 [0, 1]

    	return t;
  	}
};
```

**特殊交互说明：**与 code 2 相同。

红线是用欧式距离参数化，绿线是平均参数化。你还可以添加其他的参数化方法。

### code 4

##### work 1（估计难度 100）

求解三次样条，两段的二阶导为 0：

有解方程的方法和线性的方法，建议都尝试下

```c++
struct B3{

    vector<double > a, b, c, d; // [0, n - 1]
    vecf x, y; // [0, n]
    int p = 0;

	void getB3(){
	    ///////////////////////////////////////////////////////////////
        // y = ax^3 + bx^2 + cx + d
        // 利用 x y 求出三次样条，把系数保存到 this -> a,b,c,d 中

    	int n = x.size() - 1;
    	assert(x.size() == y.size());
    	a.resize(n);
        b.resize(n);
        c.resize(n);
        d.resize(n);
    	
        ......
}
```

### code 5

##### work 1（估计难度 10）

计算 Bernstein :

```c++
// 补全 Bernstein 中的 get() 函数
struct Bernstein{
	static double get(int n, int k, double t){
    	assert(0 <= k && k <= n);
    	/////////////////////////////////////
    	// 用定义计算 Bernstein
        
  	}
};
```

##### work 5（估计难度 10）

实现 De-Casteljau 算法：

```c++
// 补全 DeCasteljau 中的 get() 函数
struct DeCasteljau{
  	static double get(vecf x, double t){
        ///////////////////////////////////////////
        // De Casteljau 算法，我这里是递归的写法

  	}
};
```

```c++
// 补全 Bezier 中的 makeBezier 函数，计算 t 时的 x 坐标，push_back 进 ret里
// 使用 DeCasteljau 或者用 Bernstein
struct Bezier{
	static vecf makeBezier(const vecf& x, int NUM){
    	int n = x.size() - 1;
    	vecf ret;
    	double dt = 1.0 / NUM;

    	for(double t = 0; t < 1 + eps; t += dt){
    		/////////////////////////////////////
    		// 计算每个点的位置
     		
            ret.push_back(......)
    	}
    	return ret;
  	}
};
```

### code 6

##### work 1（估计难度 90）

实现 B 样条

$P(t) = \sum N_K^j (t) p_j $

```c++
// 补全 B_spline 中的 getNvec 函数
struct B_spline{
  	static vecf getNvec(int K, double t, int m){
        ////////////////////////////////////////////
        // 获取 B 样条系数 
        // m 是 平面上点的个数，K是阶数
        // 返回一个向量，表示 t 时，平面上每个点对应的 N
    }
}
```

蓝色的线是二阶B样条，两个点控制一段曲线，即折线。

红色的线是三阶B样条，三个点控制一段曲线。

黄色的线是四阶B样条，四个点控制一段曲线。

**特殊交互说明：**鼠标右键点击某个点会选中这个点，变成绿色，此时用键盘上下键可以增加或减少该点的重数。

### code 7

##### work 1（估计难度 10）

实现两种细分曲线的生成方法

- 逼近型细分：Chaiukin 方法（二次 B 样条，红色），三次 B 样条细分方法（绿色）
- 插值型细分：4 点细分方法（黄色）

```c++
// 补全 CurveSubdivision 里的函数
struct CurveSubdivision{
  	static vecf make2B(const vecf& v, int t){
        if(t <= 0){
          	return v;
        }
        //////////////////////////////////
        // Chaiukin 方法（二次 B 样条）逼近
        // v 是顶点坐标，需要细分 t 次
    
  	}

  	static vecf make3B(const vecf& v, int t){
        if(t <= 0){
          	return v;
        }
        //////////////////////////////////
        // 三次 B 样条逼近
        // v 是顶点坐标，需要细分 t 次
  	}

  	static vecf make4F(const vecf& v, int t){
        if(t <= 0){
          	return v;
        }
        //////////////////////////////////
        // 4 点细分方法
        // v 是顶点坐标，需要细分 t 次

  	}
};
```

**特殊交互说明：**用键盘上下键可以增加或减少细分的次数。

### END

也许会有后边的作业内容……