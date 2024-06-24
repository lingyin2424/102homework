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



double gauss(double x, double mu, double sigma) {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma));
}

#define circleR 7
bool same(int x1, int y1, int x2, int y2, int r = circleR){
  return std::abs(x1 - x2) <= r && std::abs(y1 - y2) <= r;
}

struct CurveSubdivision{
  static vecf make2B(const vecf& v, int t){
    if(t <= 0){
      return v;
    }
//////////////////////////////////
// Chaiukin 方法（二次 B 样条）逼近
// v 是顶点坐标，需要细分 t 次
    
    int n = v.size();
    vecf ret(n * 2);

    for(int i = 0; i < n; i ++){
      ret[2 * i] = (v[(i - 1 + n) % n] + 3 * v[i]) / 4;
      ret[2 * i + 1] = (3 * v[i] + v[(i + 1) % n]) / 4;
    }

    if(t == 1){
      return ret;
    }
    else{
      return make2B(ret, t - 1);
    }
  }

  static vecf make3B(const vecf& v, int t){
    if(t <= 0){
      return v;
    }
//////////////////////////////////
// 三次 B 样条逼近
// v 是顶点坐标，需要细分 t 次
    int n = v.size();
    vecf ret(n * 2);

    for(int i = 0; i < n; i ++){
      ret[2 * i] = (v[(i - 1 + n) % n] + 6 * v[i] + v[(i + 1) % n]) / 8;
      ret[2 * i + 1] = (4 * v[i] + 4 * v[(i + 1) % n]) / 8;
    }

    if(t == 1){
      return ret;
    }
    else{
      return make3B(ret, t - 1);
    }
  }

  static vecf make4F(const vecf& v, int t){
    if(t <= 0){
      return v;
    }
//////////////////////////////////
// 4 点细分方法
// v 是顶点坐标，需要细分 t 次

    const double alpha = 0.125;
    int n = v.size();
    vecf ret(n * 2);

    for(int i = 0; i < n; i ++){
      ret[2 * i] = v[i];
      ret[2 * i + 1] = (v[i] + v[(i + 1) % n]) / 2 + 
        alpha * ((v[i] + v[(i + 1) % n]) / 2 - (v[(i - 1 + n) % n] + v[(i + 2) % n]) / 2);
    }

    if(t == 1){
      return ret;
    }
    else{
      return make4F(ret, t - 1);
    }
  }
};

struct Poly{
  int n = 0;
  vector<float > x, y;

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
};

class MapWidget : public QWidget{
public:

  int subTime = 1;
  Poly poly;
  std::queue<std::shared_ptr<Figure > > onScreen, waitPaint;

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
    update();
    qApp->processEvents();
  }

protected:
  QColor toQColor(Vector4i color){
    return QColor(color.x(), color.y(), color.z());
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

  void drawPoly(const Poly& poly){
    if(poly.n >= 2){
      const int NUM = 3e2;
      auto x = CurveSubdivision::make2B(poly.x, subTime);
      auto y = CurveSubdivision::make2B(poly.y, subTime);

      assert(x.size() == y.size());
      int n = x.size();
      for(int i = 0; i < n; i ++){
        waitPaint.push(std::make_shared<Line >(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n], red));
      }
    }

    if(poly.n >= 2){
      auto x = CurveSubdivision::make3B(poly.x, subTime);
      auto y = CurveSubdivision::make3B(poly.y, subTime);
      
      assert(x.size() == y.size());
      int n = x.size();
      for(int i = 0; i < n; i ++){
        waitPaint.push(std::make_shared<Line >(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n], green));
      }
    }

    if(poly.n >= 2){
      auto x = CurveSubdivision::make4F(poly.x, subTime);
      auto y = CurveSubdivision::make4F(poly.y, subTime);
      
      assert(x.size() == y.size());
      int n = x.size();
      for(int i = 0; i < n; i ++){
        waitPaint.push(std::make_shared<Line >(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n], yellow));
      }
    }

    for(int i = 0; i < poly.n; i ++){
      waitPaint.push(std::make_shared<Circle >(poly.x[i], poly.y[i], red, circleR));  
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
      subTime ++;
      drawPoly(poly);
      paintIm();
    }
    else if(event -> key() == Qt::Key_Down){
      if(subTime){
        subTime --;
        drawPoly(poly);
        paintIm();
      }
    }
  }

private:
  static constexpr int pixelSize = 3; // 像素大小
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
  QApplication app(argc, argv);
  MapWidget mapWidget;
  const int WW = 2000, HH = 1000;
  mapWidget.resize(WW, HH);
  mapWidget.show();

  return app.exec();
}





