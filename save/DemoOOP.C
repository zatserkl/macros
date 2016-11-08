#include <TROOT.h>
#include <TArc.h>
#include <TMarker.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>

using std::cout;		using std::endl;

// Demo on three basic features of the Object Oriented Programming:
// -- encapsulation
// -- inheritance
// -- polymorphism

void CanvasModifiedUpdate()
{
   if (gROOT->GetListOfCanvases()->GetEntries() == 0) {
      cout<< "***Error CanvasModifiedUpdate: no TCanvas object available" <<endl;
      return;
   }
   gPad->Modified();
   gPad->Update();
}

class Location {
protected:
   Float_t x_;
   Float_t y_;
public:
   Location(Float_t x, Float_t y): x_(x), y_(y) {
      cout<< "Location::Location" <<endl;
   }
   void Print() const {
      cout<< "Location: x_ = " << x_ << " y_ = " << y_ <<endl;
   }
};

class Point: public Location {
protected:
   static const Int_t color_hide_ = 0;		// NB: assign value to static const field
private:
   TMarker* marker_;
   Int_t marker_style_;
protected:
   Int_t color_show_;
   Int_t color_curr_;
public:
   Point(Float_t x, Float_t y, Int_t color): Location(x,y)
   					     , marker_style_(20)
   					     , color_show_(color)
   					     , color_curr_(color_hide_)
   {
      cout<< "Point::Point" <<endl;
      marker_ = new TMarker(x_, y_, marker_style_);
      marker_->SetMarkerColor(color_hide_);
   }
   virtual ~Point() {
      cout<< "Point::~Point" <<endl;
      delete marker_;
   }
   void Print() const {
      cout<< "Point: x_ = " << x_
	 << " y_ = " << y_
	 << " color_show = " << color_show_
	 << " IsVisible = " << IsVisible()
	 <<endl;
   }
   // void Show()
   virtual void Show()				// polymorphism: try to remove qualifier 'virtual'
   {
      color_curr_ = color_show_;
      marker_->SetMarkerColor(color_curr_);
      marker_->DrawMarker(x_,y_);
      CanvasModifiedUpdate();
   }
   // void Hide()
   virtual void Hide()
   {
      color_curr_ = color_hide_;
      marker_->SetMarkerColor(color_curr_);
      marker_->DrawMarker(x_,y_);
      CanvasModifiedUpdate();
   }
   Bool_t IsVisible() const {return color_curr_ != color_hide_;}
   void MoveTo(Float_t x, Float_t y) {
      // method MoveTo uses methods Hide() and Show() to move the object.
      // If methods Hide() and Show() was NOT declared virtual,
      // the MoveTo will hide a Point in the old location and draw a Point in the new one.
      // If methods Hide() and Show() decalared virtual,
      // the MoveTo will use Hide() and Show() from the actual object (e.g. Circle)
      Hide();
      x_ = x;
      y_ = y;
      Show();
   }
   void Blink(Int_t nrepeat=10, Float_t delay_ms=500) {
      for (int i=0; i<nrepeat; i++) {
         Hide();
         gSystem->Sleep(delay_ms);
         Show();
         gSystem->Sleep(delay_ms);
      }
   }
};

class Circle: public Point {
protected:
   TArc* arc_;
   static const Int_t arc_line_width_ = 3;
protected:
   Float_t r_;
public:
   Circle(Float_t x, Float_t y, Float_t r, Int_t color): Point(x,y,color), r_(r)
   {
      cout<< "Circle::Circle" <<endl;
      arc_ = new TArc(x_, y_, r_);
      arc_->SetLineWidth(arc_line_width_);
      arc_->SetFillStyle(0);
      arc_->SetLineColor(color_hide_);
   }
   virtual ~Circle() {
      cout<< "Circle::~Circle" <<endl;
      delete arc_;
   }
   void Print() const {
      cout<< "Circle: x_ = " << x_
	 << " y_ = " << y_
	 << " r_ = " << r_
	 << " color_show_ = " << color_show_
	 << " IsVisible = " << IsVisible()
	 <<endl;
   }
   void Show() {							// Show inherits qualifier virtual from Point::Show
      color_curr_ = color_show_;
      arc_->SetLineColor(color_curr_);
      arc_->DrawArc(x_,y_,r_, 0, 360);
      CanvasModifiedUpdate();
   }
   void Hide() {							// Hide inherits qualifier virtual from Point::Hide
      color_curr_ = color_hide_;
      arc_->SetLineColor(color_curr_);
      arc_->DrawArc(x_,y_,r_, 0, 360);
      CanvasModifiedUpdate();
   }
};

void DemoOOP()
{
   // create a canvas

   new TCanvas("demo", "demo", 2);
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->DrawFrame(-250,-250, 250,250);

   // create (instantiate) the objects

   Point* shapes[100];     // array of pointers to Point (base object)
   Int_t nshapes = 0;

   cout<< "\nCreate a Point" <<endl;

   Point point(0,0,kRed);
   point.Print();
   // point.Show();
   shapes[nshapes++] = &point;

   cout<< "\nCreate a Circle" <<endl;

   Circle circle(100,100,50,kBlue);
   circle.Print();
   // circle.Show();
   shapes[nshapes++] = &circle;     // assign to pointer to the parent class Point

   // polymorphism: plot all the objects calling virtual method Show
   for (int ishape=0; ishape<nshapes; ishape++) {
      shapes[ishape]->Show();
   }

   // move objects

   for (int ishape=0; ishape<nshapes; ishape++)
   {
      Point* shape = shapes[ishape];
      cout<< "\nMove blinking object: Enter new x and y coordinates (space separated)" <<endl;
      shape->Blink();
      Float_t x, y;
      std::cin >> x >> y;
      shape->MoveTo(x,y);
      shape->Print();
   }

   cout<< "\nGoodbye!" <<endl;
}
