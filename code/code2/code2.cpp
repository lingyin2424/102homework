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

double gauss(double x, double mu, double sigma) {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma));
}

struct Poly{
  vector<double > x, y;
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
////////////////////////////////////////////////////// 
// 用 this -> x 和 this -> y 拟合一个 m 次多项式
// lamda 是岭回归系数，你可以先假设 lamda = 0，即没有岭回归
// 返回拟合的多项式的系数，小的下标对应低次的系数
    to do
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
  double lamdaNum = 1e-2;

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
    cout << "fitNum = " << fitNum << " lamdaNum = " << lamdaNum << endl;

    if(poly.x.size() >= 2){
      {
        int m = fitNum;
        auto a = poly.polynomialFit(m);
        int lastx = 0;
        int lasty = poly.getPoly(a, 0);

        for(int x = 1; x < size().width(); x += pixelSize){
          int y = poly.getPoly(a, x);
          waitPaint.push(std::make_shared<Line >(lastx, lasty, x, y, green));
          lastx = x;
          lasty = y;
        }
      }

      {
        int m = fitNum;
        auto a = poly.polynomialFit(m, lamdaNum);
        int lastx = 0;
        int lasty = poly.getPoly(a, 0);

        for(int x = 1; x < size().width(); x += pixelSize){
          int y = poly.getPoly(a, x);
          waitPaint.push(std::make_shared<Line >(lastx, lasty, x, y, yellow));
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






