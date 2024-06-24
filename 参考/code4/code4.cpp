#include <iostream>
#include <QtWidgets>
#include <Eigen/Dense>
#include <vector>
#include <assert.h>
#include <utility>
#include <memory>
#include <tuple>
#include <queue>

const double eps = 1e-3;

using veci = std::vector<int >;
using vecf = std::vector<float >;

template<class T >
T pp(T x){
    return x * x;
}

template<class T >
T ppp(T x){
    return x * x * x;
}

using Eigen::Matrix4f;
using Eigen::Matrix4i;
using Eigen::Vector3f;
using Eigen::Vector3i;
using Eigen::Vector4f;
using Eigen::Vector4i;
using std::vector;

using std::cout;
using std::endl;

using Image = vector<vector<Vector4i >>;
using std::max;
using std::min;
using std::pair;
const Vector4i black(0, 0, 0, 0);
const Vector4i white(255, 255, 255, 0);
const Vector4i red(255, 0, 0, 0);
const Vector4i orange(255, 128, 128, 0);
const Vector4i blue(0, 0, 255, 0);
const Vector4i green(30, 255, 0, 0);
const Vector4i yellow(255, 255, 0, 0);


class Figure{ 
  
  virtual std::tuple<veci, veci, vector<Vector4i> > makePaintList(int ww, int hh, int pixelSize) = 0;
public:
  virtual ~Figure(){};  // 虚析构函数

  std::tuple<veci, veci, vector<Vector4i> > getPaintList(int ww, int hh, int pixelSize){
    auto [x, y, c] = makePaintList(ww, hh, pixelSize);
    assert(x.size() == y.size() && y.size() == c.size());

    for(int & _x : x){
      _x = max(0, min(ww - 1, _x));
    }
    for(int & _y : y){
      _y = max(0, min(hh - 1, _y));
    }

    return {x, y, c};
  }
};

class Line : public Figure{
  int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
  Vector4i color = white;
public:
  Line(int _x0, int _y0, int _x1, int _y1, Vector4i _color)
    : x0{_x0}, y0{_y0}, x1{_x1}, y1{_y1}, color{_color} {}

  std::tuple<veci, veci, vector<Vector4i> > makePaintList(int ww, int hh, int pixelSize) override {  
    veci retx, rety;  
    vector<Vector4i > retColor;
    using std::swap;
    if(abs(x0 - x1) >= abs(y0 - y1)){
      if(x0 > x1){
        swap(x0, x1);
        swap(y0, y1);
      }
      int len = x1 - x0 + 1;
      for(int x = x0; x <= x1; x += pixelSize){
        int y = (y0 * (x1 - x) + y1 * (x - x0 + 1)) / len;
        retx.push_back(x);
        rety.push_back(y);
        retColor.push_back(color);
      }
    }
    else{
      if(y0 > y1){
        swap(x0, x1);
        swap(y0, y1);
      }
      int len = y1 - y0 + 1;
      for(int y = y0; y <= y1; y += pixelSize){
        int x = (x0 * (y1 - y) + x1 * (y - y0 + 1)) / len;
        retx.push_back(x);
        rety.push_back(y);
        retColor.push_back(color);
      }
    }

    return {retx, rety, retColor};
  }

  ~Line() override {};
};


class Circle : public Figure{
  int x = 0, y = 0, r0 = 0, r = 0;
  Vector4i color = white;
public:
  Circle(int _x, int _y, Vector4i _color, int _r, int _r0 = 0)
    : x{_x}, y{_y}, color{_color}, r{_r}, r0{_r0} {}

  std::tuple<veci, veci, vector<Vector4i> > makePaintList(int ww, int hh, int pixelSize) override { 
    veci retx, rety;  
    vector<Vector4i > retColor;

    for(int i = x - r / pixelSize * pixelSize; i <= x + r; i += pixelSize){
      for(int j = y - r / pixelSize * pixelSize; j <= y + r; j += pixelSize){
        if(pp(i - x) + pp(j - y) <= pp(r)
        && pp(i - x) + pp(j - y) >= pp(r0)){
          retx.push_back(i);
          rety.push_back(j);
          retColor.push_back(color);
          
        }
      }
    }

    return {retx, rety, retColor};
  }
  ~Circle() override {};
};

#define circleR 7
bool same(int x1, int y1, int x2, int y2, int r = circleR){
  return std::abs(x1 - x2) <= r && std::abs(y1 - y2) <= r;
}

struct Parameterization{
  // 对曲线的时间做参数化

  static vector<float > distance(const vecf& x, const vecf& y){

    int n = x.size();
    assert(n >= 2);
    assert(x.size() == y.size());
    vecf t(n);
    
    for(int i = 1; i < n; i ++){
      t[i] = t[i - 1] + sqrt(pp(x[i] - x[i - 1]) + pp(y[i] - y[i - 1]));
    }
    for(int i = 1; i < n; i ++){
      t[i] /= t.back();
    }

    return t;
  }

  static vector<float > average(const vecf& x, const vecf& y){

    int n = x.size();
    assert(n >= 2);
    assert(x.size() == y.size());
    vecf t(n);
    for(int i = 0; i < n; i ++){
      t[i] = (double)i / (n - 1);
    }

    return t;
  }
};
double gauss(double x, double mu, double sigma) {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma));
}

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
    

    /* 这里是用矩阵求解的方法，建议先尝试
    cout << " n = " << n << endl;
    Eigen::MatrixXd A(4 * n, 4 * n);
    Eigen::VectorXd Y(4 * n);

    for(int i = 0; i < n * 4; i ++){
      for(int j = 0; j < n * 4; j ++){
        A(i, j) = 0;
      }
    }

    for(int i = 0, j = 0; i < n; i ++){
      A(i, j ++) = ppp(x[i]); 
      A(i, j ++) = pp(x[i]); 
      A(i, j ++) = x[i];
      A(i, j ++) = 1;
    }

    for(int i = 0, j = 0; i < n; i ++){
      A(i + n, j ++) = ppp(x[i + 1]); 
      A(i + n, j ++) = pp(x[i + 1]); 
      A(i + n, j ++) = x[i + 1];
      A(i + n, j ++) = 1;
    }

    {// S0''(x0) == 0
     // 6a_0x_0 + 2b_0 == 0
      A(n + n, 0) = 6 * x[0];
      A(n + n, 1) = 2;

      // S_{n - 1}''(x_n) == 0
      // 6 a_{n - 1}x_{n - 1} + 2 b_{n - 1} == 0

      A(n + n + n, (n - 1) << 2) = 6 * x[n];
      A(n + n + n, (n - 1) << 2 | 1) = 2;      
    }

    for(int i = 1, j1 = 0, j2 = 4; i < n; i ++){
      int id = i + n + n;
      A(id, j1 ++) = -3 * pp(x[i]);
      A(id, j1 ++) = -2 * x[i];
      A(id, j1 ++) = -1;
      A(id, j1 ++) = 0;

      A(id, j2 ++) = +3 * pp(x[i]);
      A(id, j2 ++) = +2 * x[i];
      A(id, j2 ++) = +1;
      A(id, j2 ++) = 0;
    }

    for(int i = 1, j1 = 0, j2 = 4; i < n; i ++){
      int id = i + n + n + n;
      A(id, j1 ++) = -6 * x[i];
      A(id, j1 ++) = -2;
      A(id, j1 ++) = 0;
      A(id, j1 ++) = 0;

      A(id, j2 ++) = +6 * x[i];
      A(id, j2 ++) = +2;
      A(id, j2 ++) = 0;
      A(id, j2 ++) = 0;
    }


    for(int i = 0; i < n; i ++){
      Y[i] = y[i];
      Y[i + n] = y[i + 1];
      Y[i + n + n] = 0;
      Y[i + n + n + n] = 0;
    }
    
    a.resize(n);
    b.resize(n);
    c.resize(n);
    d.resize(n);

    Eigen::VectorXd B = A.fullPivHouseholderQr().solve(Y);
    for(int i = 0; i < n; i ++){
      a[i] = B[i << 2];
      b[i] = B[i << 2 | 1];
      c[i] = B[i << 2 | 2];
      d[i] = B[i << 2 | 3];
    }
    */


    // 这里是线性的方法
    vector<double > cc(n), cp(n), dd(n), dp(n);

    for(int i = 0; i < n; i ++){
      cc[i]  = - pp(x[i + 1] - x[i]) / 6;
      cp[i] = y[i + 1];
      dd[i]  = - pp(x[i + 1] - x[i]) / 6;
      dp[i] = y[i    ];
    }

    Eigen::VectorXd m = Eigen::VectorXd::Zero(n + 1);
    m[0] = m[n] = 0;

    if(n > 1){
      Eigen::VectorXd Y(n - 1);

      for(int i = 0; i <= n - 2; i ++){
        auto E = x[i + 2] - x[i + 1];
        auto F = x[i + 1] - x[i];
        Y[i] = (-cp[i] * E + dp[i] * E + cp[i + 1] * F - dp[i + 1] * F);
      }

      vector<double > _a(n - 1), _b(n - 1), _c(n - 1);
      for(int i = 0; i <= n - 2; i ++){
        auto E = x[i + 2] - x[i + 1];
        auto F = x[i + 1] - x[i];

        if(i != 0){
          _a[i] = -dd[i] * E;
        }

        _b[i] = pp(x[i + 1] - x[i    ]) / 2.0 * E 
              + pp(x[i + 2] - x[i + 1]) / 2.0 * F
              + (cc[i] * E + dd[i + 1] * F);

        if(i != n - 2){
          _c[i] = -cc[i + 1] * F;
        }
      }

      for(int i = 1; i <= n - 2; i ++){
        auto k = _a[i] / _b[i - 1];
        _a[i] -= k * _b[i - 1];
        _b[i] -= k * _c[i - 1];
         Y[i] -= k *  Y[i - 1];
      }

      for(int i = 0; i <= n - 2; i ++){
        assert(std::abs(_a[i]) < eps);
      }
      
      for(int i = n - 3; i >= 0; i --){
        auto k = _c[i] / _b[i + 1];
        _c[i] -= k * _b[i + 1];
         Y[i] -= k *  Y[i + 1];
      }

      for(int i = 0; i <= n - 2; i ++){
        assert(std::abs(_c[i]) < eps);
      }

      for(int i = 0; i <= n - 2; i ++){
        m[i + 1] = Y[i] / _b[i];
      }
    }

    Eigen::VectorXd C(n), D(n);

    a.resize(n);
    b.resize(n);
    c.resize(n);
    d.resize(n);
 
    for(int i = 0; i < n; i ++){
      C[i] = cc[i] * m[i + 1] + cp[i];
      D[i] = dd[i] * m[i    ] + dp[i];
    
      a[i] = (m[i + 1] / 6 * (+1)             + m[i] / 6 * (-1)                );
      b[i] = (m[i + 1] / 6 * (-3 * x[i])      + m[i] / 6 * (+3 * x[i + 1])     );
      c[i] = (m[i + 1] / 6 * (+3 * pp(x[i]))  + m[i] / 6 * (-3 * pp(x[i + 1])) ) + C[i] - D[i];
      d[i] = (m[i + 1] / 6 * (-1 * ppp(x[i])) + m[i] / 6 * (+1 * ppp(x[i + 1])))            
           - x[i] * C[i] + x[i + 1] * D[i];


      a[i] /= (x[i + 1] - x[i]);
      b[i] /= (x[i + 1] - x[i]);
      c[i] /= (x[i + 1] - x[i]);
      d[i] /= (x[i + 1] - x[i]);
    }


  }

  double get(double t){
    while(p != x.size() - 2 && t > x[p + 1]){
      p ++;
    }
    return a[p] * ppp(t) + b[p] * pp(t) + c[p] * t + d[p];
  }
};


struct Poly{
  vecf x, y;
  int n = 0;

  void add(double _x, double _y){
    x.push_back(_x);
    y.push_back(_y);
    n ++;
  }
  
  void remove(int id){
    assert(id < n);
    x.erase(x.begin() + id);
    y.erase(y.begin() + id);
    n --;
  }

  int find(int _x, int _y){
    assert(n == x.size());
    for(int i = 0; i < n; i ++){
      if(same(_x, _y, x[i], y[i])){
        return i;
      }
    }

    return -1;
  }

  vector<double> polynomialFit(int m, double lamda = 0){
    int n = x.size();
    Eigen::MatrixXd A(m, m);
    for(int i = 0; i < m; i ++){
      for(int j = 0; j < m; j ++){
        double tmp = lamda;
        for(int v = 0; v < n; v ++){
          tmp += pow(x[v], i) * pow(x[v], j);
        }
        A(i, j) = tmp;
      }
    }

    Eigen::VectorXd Y(m);
    for(int i = 0; i < m; i ++){
      double tmp = 0;
      for(int v = 0; v < n; v ++){
        tmp += pow(x[v], i) * y[v];
      }
      Y[i] = tmp;
    }

    Eigen::VectorXd ret = A.partialPivLu().solve(Y);
    return vector<double> (ret.data(), ret.data() + ret.size());
  }

  double getPoly(const vector<double > a, double x){
    int n = a.size();
    double ret = 0;
    for(int i = n - 1; i >= 0; i --){
      ret = (ret * x + a[i]);
    }
    return ret;
  }
};

class MapWidget : public QWidget{
public:
  Poly poly;
  std::queue<std::shared_ptr<Figure > > onScreen, waitPaint;
  int fitNum = 3; 
  double lamdaNum = 1e-3;

  struct{
    int x1, y1, x2, y2;
    void setIn(int _x, int _y){
      x1 = _x;
      y1 = _y;
    }
    
    void setOut(int _x, int _y){
      x2 = _x;
      y2 = _y;
    } 

    void update(Poly& poly){
      int idx1 = poly.find(x1, y1);
      int idx2 = poly.find(x2, y2);
      
      if(idx1 != -1){
        if(idx2 == -1){
          poly.x[idx1] = x2;
          poly.y[idx1] = y2;
        }
        else if(idx1 == idx2){
          poly.remove(idx1);
        }
        else{}
      }
      else{
        if(idx2 == -1){
          poly.add(x2, y2);
        }
        else{}
      }
    }
  }manageClick;

  void paintIm(){
    // 立即重绘
    update();
    qApp->processEvents();
  }

protected:
  QColor toQColor(Vector4i color){
    return QColor(color.x(), color.y(), color.z());
  }


  void drawPoly(Poly& poly){
    if(poly.x.size() >= 2){
      
      {
        auto t = Parameterization::distance(poly.x, poly.y);
        Poly polyx, polyy;

        B3 b3x, b3y;
        b3x.x = t;
        b3x.y = poly.x;
        b3x.getB3();

        b3y.x = t;
        b3y.y = poly.y;
        b3y.getB3();


        int lastx = b3x.get(0);
        int lasty = b3y.get(0);
        for(float i = 0; i < 1 + eps; i += eps){
          int x = b3x.get(i);
          int y = b3y.get(i);
          waitPaint.push(std::make_shared<Line >(lastx, lasty, x, y, red));
          lastx = x;
          lasty = y;
        }
      }
    }

    for(int i = 0; i < poly.n; i ++){
      waitPaint.push(std::make_shared<Circle >(poly.x[i], poly.y[i], red, circleR));  
    }
  }

  void mousePressEvent(QMouseEvent *event) override {
    if(event->button() == Qt::RightButton){}
    else if(event->button() == Qt::LeftButton){
      QPoint clickPosition = event->pos();  
      int clickX = clickPosition.x();
      int clickY = clickPosition.y();
      manageClick.setIn(clickX, clickY);
    }
  }

  void mouseReleaseEvent(QMouseEvent *event) override {
    if(event->button() == Qt::RightButton){}
    else if(event->button() == Qt::LeftButton){
      QPoint clickPosition = event->pos();  
      int clickX = clickPosition.x();
      int clickY = clickPosition.y();

      manageClick.setOut(clickX, clickY);
      manageClick.update(poly);
      drawPoly(poly);
      paintIm();
    }
  }

  void mouseMoveEvent(QMouseEvent * event) override{
      QPoint clickPosition = event->pos();  
      int clickX = clickPosition.x();
      int clickY = clickPosition.y();
      int t1 = poly.find(manageClick.x1, manageClick.y1);  
      
      if(t1 != -1){
        auto tmpPoly = poly;
        tmpPoly.x[t1] = clickX;
        tmpPoly.y[t1] = clickY;
        drawPoly(tmpPoly);
        paintIm();
      }
  }

  void paintEvent(QPaintEvent *event) override{
    Q_UNUSED(event);
    QPainter painter(this);
    int ww = size().width();
    int hh = size().height();

    while(waitPaint.empty() == false){
      auto [x, y, c] = waitPaint.front() -> getPaintList(ww, hh, pixelSize);
      int n = x.size();
      for(int i = 0; i < n; i ++){
        painter.setBrush(QColor(c[i].x(), c[i].y(), c[i].z()));
        painter.drawRect(x[i], y[i], pixelSize, pixelSize);
      }

      onScreen.push(waitPaint.front());
      waitPaint.pop();
    }
  }

  void keyPressEvent(QKeyEvent *event) override {
    if(event -> key() == Qt::Key_Up){
      lamdaNum *= 2;
      drawPoly(poly);
      paintIm();
    }
    else if(event -> key() == Qt::Key_Down){
      lamdaNum /= 2;
      drawPoly(poly);
      paintIm();
    }
    else if(event -> key() == Qt::Key_Left){
      if(fitNum > 2){
        fitNum --;
      }
      drawPoly(poly);
      paintIm();
    }
    else if(event -> key() == Qt::Key_Right){
      fitNum ++;
      drawPoly(poly);
      paintIm();
    }
  }

private:
  static constexpr int pixelSize = 3; // 像素大小
};

int getRand(int l, int r){
  return rand() % (r - l + 1) + l;
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
  QApplication app(argc, argv);
  MapWidget mapWidget;
  const int WW = 2000, HH = 1000;
  mapWidget.resize(WW, HH);
  mapWidget.show();

  return app.exec();
}
