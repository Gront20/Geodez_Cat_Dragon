//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>

//---------------------------------------------------------------------------
class Complex
{
public:
    bool IsInfinity;
    double re, im;
    Complex();
    Complex(double re, double im=0);
    Complex(const Complex& other);
    Complex operator+(const Complex &c) const;
    Complex operator-(const Complex &c) const;
    Complex operator/(const Complex &c) const;
    Complex operator^(const int &c) const;
    Complex operator*(const Complex& other) const;
    Complex operator*(const double& other) const;
    Complex& operator*=(const Complex& other);
    Complex& operator=(const Complex& other);
    bool operator== (const Complex& other) const;
    Complex operator~() const;
    Complex operator-() const;
    double Arg();
    double Abs();
    Complex sqrt();
    Complex revers();
    Complex conj();
    Complex Polar(double mod, double arg);
};

typedef struct {
  TList *polygon;
  TList *data;
  TList *grid;
  Complex d0;
  Complex d1;
  Complex od0;
  Complex ocenter;
  Complex center;
} ConformalMapRec;


typedef struct {
  TList *polygon;
  TList *data;
  TList *grid;
  Complex d0;
  Complex d1;
  Complex od0;
  Complex ocenter;
  Complex center;
} InverseConformalMapRec;

//---------------------------------------------------------------------------

class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TButton *Button1;
        TButton *Button2;
        TButton *Button3;
        TButton *Button4;
        TPageControl *PageControl1;
        TTabSheet *TabSheet1;
        TTabSheet *TabSheet2;
        TTabSheet *TabSheet3;
        TTabSheet *TabSheet4;
        TTabSheet *TabSheet5;
        TPanel *Panel1;
        TImage *Image1;
        TImage *Image2;
        TImage *Image3;
        TImage *Image4;
        TImage *Image5;
        TButton *Button5;
        TButton *Button6;
        TTabSheet *TabSheet6;
        TButton *Button7;
        TImage *Image6;
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall Button3Click(TObject *Sender);
        void __fastcall Button4Click(TObject *Sender);
        void __fastcall Button5Click(TObject *Sender);
        void __fastcall Button6Click(TObject *Sender);
        void __fastcall Button7Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
        TList *test;
        ConformalMapRec *conmap,*conmap1;
        InverseConformalMapRec *invconmap,*invconmap1;
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
 