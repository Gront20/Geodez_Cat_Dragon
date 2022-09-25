//---------------------------------------------------------------------------
 
#include <vcl.h>
#pragma hdrstop
#include <math.h>
#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

#define COUNT_PER_SHAPE 100

//---------------------------------------------------------------------------
double __fastcall MaxValuePoint(TList *points)
{
   double m=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      if(fabs(t1->re)>m) m=fabs(t1->re);
      if(fabs(t1->im)>m) m=fabs(t1->im);
   }
   return m;
}

//---------------------------------------------------------------------------
double __fastcall MaxLine(TList *points)
{
   double m=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      int l=i+1;
      if(l>d-1) l=d-i-1;
      Complex *t2=(Complex *)(points->Items[l]);
      Complex t3=*t1-*t2;
      if(m<t3.Abs()) m=t3.Abs();
   }
   return m;
}

//---------------------------------------------------------------------------
double __fastcall Razmer(TList *points)
{
   double m=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      for(int j=0;j<d;j++)
      {
         Complex *t2=(Complex *)(points->Items[j]);
         double v=sqrt(pow((t1->re-t2->re),2)+pow((t1->im-t2->im),2));
         if(m<v) m=v;
      }
   }
   return m;
}

//---------------------------------------------------------------------------
bool new_inside(double x,double y,TList *polygon)
{
    static const int q_patt[2][2]= { {0,1}, {3,2} };
    if (polygon->Count<3) return false;

    Complex pred_pt=*(Complex *)(polygon->Items[0]);
    pred_pt.re-=x;
    pred_pt.im-=y;

    int pred_q=q_patt[pred_pt.im<0][pred_pt.re<0];

    int w=0;

    for(int i=0;i<polygon->Count;i++)
    {
       Complex *temp=(Complex *)(polygon->Items[i]);
       Complex cur_pt=*temp;
       cur_pt.re-=x;
       cur_pt.im-=y;

       int q=q_patt[cur_pt.im<0][cur_pt.re<0];
       switch (q-pred_q)
       {
           case -3:++w;break;
           case 3:--w;break;
           case -2:if(pred_pt.re*cur_pt.im>=pred_pt.im*cur_pt.re) ++w;break;
           case 2:if(!(pred_pt.re*cur_pt.im>=pred_pt.im*cur_pt.re)) --w;break;
       }
       pred_pt = *temp;
       pred_q = q;
    }
    return w!=0;
}

//---------------------------------------------------------------------------

int inside(double px,double py,TList *lst) {
    int i,j,s;
    int n=lst->Count;
    s=0; j=n-1;
    for(i=0;i<n;j=i++) {
        Complex *t1=(Complex *)(lst->Items[i]);
        double xi=t1->re;
        double yi=t1->im;
        t1=(Complex *)(lst->Items[j]);
        double xj=t1->re;
        double yj=t1->im;
        if ( (py<yi ^ py<yj) || py==yi || py==yj ) {
            if ( (px<xi ^ px<xj) || px==xi || px==xj ) {
                if ( yi<yj ^ (xi-px)*(yj-py)>(yi-py)*(xj-px) ) s^=1;
            } else {
                if (px>xi && yi!=yj) s^=1;
            }
        }
    }
    if(s) s=new_inside(px,py,lst);
    return s;
}

//---------------------------------------------------------------------------
int __fastcall sign(double val)
{
   if(val<0) return -1;
   else return 1;
}

//---------------------------------------------------------------------------
Complex::Complex()
{
   IsInfinity=false;
   re=0.0;
   im=0.0;
}

//---------------------------------------------------------------------------
Complex::Complex(double fre, double fim)
{
   IsInfinity=false;
   re=fre;
   im=fim;

}
//---------------------------------------------------------------------------
Complex::Complex(const Complex& x)
{
   IsInfinity=x.IsInfinity;
   re=x.re;
   im=x.im;
}
//---------------------------------------------------------------------------
Complex Complex::operator*(const Complex& other) const
{
    Complex temp;
    temp.re=re*other.re-im*other.im;
    temp.im=re*other.im+im*other.re;
    return temp;
}

//---------------------------------------------------------------------------
Complex Complex::operator*(const double& other) const
{
    return Complex(re*other, im*other);
}

//---------------------------------------------------------------------------
Complex& Complex::operator*=(const Complex& other)
{
    Complex temp=*this;
    re=temp.re*other.re - temp.im*other.im;
    im=temp.re*other.im + temp.im*other.re;
    return *this;
}

//---------------------------------------------------------------------------
Complex Complex::operator-() const
{
    return Complex(-re, -im);
}

//---------------------------------------------------------------------------
Complex Complex::operator~() const
{
    return Complex(re, -im);
}

//---------------------------------------------------------------------------
bool Complex::operator== (const Complex& other) const
{
    return (re == other.re && im == other.im);
}

//---------------------------------------------------------------------------
Complex& Complex::operator=(const Complex& other)
{
    if(this != &other)
    {
        re=other.re;
        im=other.im;
        IsInfinity=other.IsInfinity;
    }
    return *this;
}

//---------------------------------------------------------------------------
double Complex::Arg()
{
   double delim=re;
   if(fabs(delim)<1e-12) delim=1e-12;
   return atan2(im,delim);
}

//---------------------------------------------------------------------------
double Complex::Abs()
{
   return ::sqrt(re * re + im * im);
}


//---------------------------------------------------------------------------
Complex Complex::operator - (const Complex &c) const
{
	  Complex temp=*this;
	  temp.re = temp.re - c.re;
	  temp.im = temp.im - c.im;

	  return temp;
}

//---------------------------------------------------------------------------
Complex Complex::operator + (const Complex &c) const
{
	  Complex temp;

	  temp.re = re + c.re;
	  temp.im = im + c.im;

	  return temp;
}
//---------------------------------------------------------------------------
Complex Complex::operator / (const Complex &c) const
{
	  Complex temp;

          double r = c.re * c.re + c.im * c.im;
          if(fabs(r)<1e-34) r=1e-34; 
          temp.re = (re * c.re + im * c.im) / r;
          temp.im = (im * c.re - re * c.im) / r;

	  return temp;
}

//---------------------------------------------------------------------------

Complex Complex::sqrt()
{
    Complex temp=*this;
    double abs=::sqrt(temp.Abs());
    return Complex(abs*cos(temp.Arg()/2),abs*sin(temp.Arg()/2));
}

//---------------------------------------------------------------------------

Complex Complex::revers()
{
    Complex temp=*this;
    double tmp;
    tmp=temp.re;
    temp.re=-temp.im;
    temp.im=tmp;
    return temp;
}

//---------------------------------------------------------------------------

Complex Complex::conj()
{
    Complex temp=*this;
    temp.im=(-1)*temp.im;
    return temp;
}


//---------------------------------------------------------------------------
Complex Complex::Polar(double mod, double arg)
{
   re=mod * cos(arg);
   im=mod * sin(arg);
   return *this;
}

//---------------------------------------------------------------------------
Complex Complex::operator^(const int &c) const
{
  Complex temp=*this;

  for(int i=0;i<c-1;i++) temp=temp*temp;
  return temp;
}
//---------------------------------------------------------------------------
Complex __fastcall f(Complex z,Complex a)
{
   Complex Result;
   if(z.IsInfinity)
   {
      Complex t1=Complex(a.re);
      Complex t2=Complex(a.im);
      Complex t3=t1/t2.revers();
      Result=((t3^2)-Complex(1)).sqrt();
      Result.IsInfinity=false;
   }
   else
   {
     if(fabs(a.im)<1e-12)
     {
         Complex t0=Complex(a.re);
         Complex t1=(z/t0)^2;
         Complex t2=Complex(1);
         Complex z1=t1-t2;
         Result=z1.sqrt();
     } else
     {
       Complex c=Complex(a.re/(a.Abs()*a.Abs()));
       Complex d=Complex(a.im/(a.Abs()*a.Abs()));
       Complex t1=(z*c)/(Complex(1)+z.revers()*d);
       Complex t2=t1^2;
       Complex t3=(t2-Complex(1)).sqrt();
       if(t3.im*t1.im<0) t3=t3*(-1.0);
       Result=t3;
     }

   }
   return Result;
}

//---------------------------------------------------------------------------
Complex __fastcall finv(Complex w,Complex a)
{
   Complex Result;
   double tol=1e-12;
   if(fabs(a.im)<tol)
   {
     Complex t2=(w^2)+Complex(1);
     Complex t3=Complex(a.re)^2;
     Complex rew=t2*t3;
     Result=rew.sqrt();
   }
   else
   {
      if(w.IsInfinity)
      {
         Complex d=Complex(a.im/(a.Abs()*a.Abs()));
         Result=Complex(-1)/d.revers();
         Result.IsInfinity=false;
      }
      else
      {
         Complex c=Complex(a.re/(a.Abs()*a.Abs()));
         Complex d=Complex(a.im/(a.Abs()*a.Abs()));
         Complex t1=((w^2)+Complex(1)).sqrt();
         Complex t2=t1/(c-t1.revers()*d);
         Result=t2;
      }
   }
   return Result;
}

//---------------------------------------------------------------------------
Complex __fastcall Composition(TList *A,int count,Complex z)
{
   Complex Result=z;
   for(int i=0;i<count;i++)
   {
      Complex *tmp=(Complex *)(A->Items[i]);
      Result=f(Result,*tmp);
   }
   return Result;
}


//---------------------------------------------------------------------------
Complex __fastcall InvComposition(TList *A,int count,Complex w)
{
   Complex Result=w;
   for(int i=0;i<count;i++)
   {
      Complex *tmp=(Complex *)(A->Items[count-i-1]);
      Result=finv(Result,*tmp);
      Result.IsInfinity=false;
   }
   return Result;

}

//---------------------------------------------------------------------------
Complex __fastcall PhiFirst(Complex z,Complex point1,Complex point2)
{
   Complex Result(0.0,0.0);
   if(z==point1) Result.IsInfinity=true;
   else
   {
      Complex t1=z-point2;
      Complex t2=z-point1;
      Complex t3=t1/t2;
      Result=t3.sqrt();
   }
   return Result;
}


//---------------------------------------------------------------------------
Complex __fastcall PhiFirstInv(Complex w,Complex point1,Complex point2)
{
   Complex Result(0.0,0.0);
   if(w.IsInfinity)
   {
      Result=point1;
      return Result;
   }
   Complex t1=w^2;
   Result=Complex((t1*point1-point2)/(t1-Complex(1.0)));
   return Result;
}


//---------------------------------------------------------------------------
Complex __fastcall PhiLast(Complex z,Complex last)
{
   Complex Result;
   double sign=-1.0;

   if(z.IsInfinity) return last^2;
   Complex t0=z/last;
   if((Complex(1.0)-t0).Abs()<1e-12)
   {
      Result.IsInfinity=true;
      return Result;
   }
   Complex t01=z/(Complex(1.0)-t0);
   Result=(t01^2)*sign;
   return Result;
}

//---------------------------------------------------------------------------
Complex __fastcall PhiLastInv(Complex w,Complex last)
{
   Complex Result(0.0,0.0);
   double sign=-1.0;

   if(w.Abs()<1e-12) return Result;
   if(w.IsInfinity) return last;
   Complex t0=(w*sign).sqrt();
   Complex t3=last+t0;
   if(t3.Abs()<1e-12)
   {
      Result.IsInfinity=true;
      return Result;
   }
   Complex t1=t0*last;
   Result=t1/t3;
   return Result;
}

//---------------------------------------------------------------------------
Complex __fastcall HalfPlaneToDisk(Complex z,Complex a)
{
   Complex t1=z-a;
   Complex t2=z-a.conj();
   Complex t3=t1/t2;
   return t3;
}


//---------------------------------------------------------------------------
Complex __fastcall DiskToHalfPlane(Complex w,Complex a)
{
   return (w*a.conj()-a)/(w-Complex(1.0));
}

//---------------------------------------------------------------------------
TList *__fastcall Zapol(TList *A)
{
    TList *dna=new TList();
    double pntcount=MaxLine(A)/COUNT_PER_SHAPE;
    int L=A->Count;
    for(int i=0;i<L;i++)
    {
      Complex *t1=(Complex *)(A->Items[i]);
      int l=i+1;
      if(l>L-1) l=L-i-1;
      Complex *t2=(Complex *)(A->Items[l]);
      Complex t3=*t1-*t2;
      int n=ceil(t3.Abs()/pntcount);
      for(int k=0;k<n;k++)
      {
        double x=t1->re+k*pntcount*(t2->re-t1->re)/t3.Abs();
        double y=t1->im+k*pntcount*(t2->im-t1->im)/t3.Abs();
        Complex *tmp=new Complex(x,y);
        dna->Add(tmp);
      }
    }
    return dna;
}

//---------------------------------------------------------------------------
void __fastcall center(TList *points, Complex &v)
{
   double rz=0;
   double iz=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      rz+=t1->re;
      iz+=t1->im;
   }
   v.re=rz/points->Count;
   v.im=iz/points->Count;
}


//---------------------------------------------------------------------------
ConformalMapRec * __fastcall GenerateConformalMap(TList *polygon,TList *data,Complex &d0,Complex &d1,Complex &od0,Complex &ocenter,Complex &center)
{
   ConformalMapRec *gnr=new ConformalMapRec;
   Complex *t1;
   gnr->polygon=new TList();
   for(int i=0;i<polygon->Count;i++)
   {
      t1=new Complex;
      t1=(Complex *)(polygon->Items[i]);
      gnr->polygon->Add(t1);
   }
   gnr->data=new TList();
   for(int i=0;i<data->Count;i++)
   {
      t1=new Complex;
      t1=(Complex *)(data->Items[i]);
      gnr->data->Add(t1);
   }
   gnr->d0=d0;
   gnr->d1=d1;
   gnr->od0=od0;
   gnr->ocenter=ocenter;
   gnr->center=center;
   return gnr;
}

//---------------------------------------------------------------------------
void __fastcall DeleteConformalMap(ConformalMapRec *val)
{
   for(int i=0;i<val->polygon->Count;i++) delete val->polygon->Items[i];
   val->polygon->Clear();
   for(int i=0;i<val->data->Count;i++) delete val->data->Items[i];
   val->data->Clear();
   delete val;
}


//---------------------------------------------------------------------------
InverseConformalMapRec * __fastcall GenerateInverseConformalMap(TList *polygon,TList *data,Complex &d0,Complex &d1,Complex &od0,Complex &ocenter,Complex center)
{
   InverseConformalMapRec *gnr=new InverseConformalMapRec;
   Complex *t1;
   gnr->polygon=new TList();
   for(int i=0;i<polygon->Count;i++)
   {
      t1=new Complex;
      t1=(Complex *)(polygon->Items[i]);
       gnr->polygon->Add(t1);
   }
   gnr->data=new TList();
   for(int i=0;i<data->Count;i++)
   {
      t1=new Complex;
      t1=(Complex *)(data->Items[i]);
      gnr->data->Add(t1);
   }
   gnr->d0=d0;
   gnr->d1=d1;
   gnr->od0=od0;
   gnr->ocenter=ocenter;
   gnr->center=center;
   gnr->center=center;
   return gnr;
}

//---------------------------------------------------------------------------
void __fastcall DeleteInverseConformalMap(InverseConformalMapRec *val)
{
   for(int i=0;i<val->polygon->Count;i++) delete val->polygon->Items[i];
   val->polygon->Clear();
   for(int i=0;i<val->data->Count;i++) delete val->data->Items[i];
   val->data->Clear();
   delete val;
}

//---------------------------------------------------------------------------
ConformalMapRec * __fastcall ConformalMap(TList *polygon,Complex center)
{
   TList *Zapolarray=Zapol(polygon);
   TList *temporary=new TList();
   Complex *d0=(Complex *)(Zapolarray->Items[0]);
   Complex *d1=(Complex *)(Zapolarray->Items[1]);
   for(int k=2;k<Zapolarray->Count;k++)
   {
      Complex *t1=new Complex();
      Complex *dk=(Complex *)(Zapolarray->Items[k]);
      *t1=Composition(temporary,temporary->Count,PhiFirst(*dk,*d0,*d1));
      temporary->Add(t1);
   }
   Complex od0=Composition(temporary,temporary->Count,PhiFirst(*d0,*d0,*d1));
   Complex ocenter=PhiLast(Composition(temporary,temporary->Count,PhiFirst(center,*d0,*d1)),od0);

   TStringList *test=new TStringList();
   for(int r=0;r<temporary->Count;r++)
   {
      Complex *rrr=(Complex *)(temporary->Items[r]);
      test->Add(FloatToStr(rrr->re)+" "+FloatToStr(rrr->im));
   }
   test->SaveToFile("out.txt");
   delete test;

   return GenerateConformalMap(Zapolarray,temporary,*d0,*d1,od0,ocenter,center);
}


//---------------------------------------------------------------------------
Complex __fastcall GetConformalMap(ConformalMapRec *CM,Complex z)
{
    TList *q=CM->data;
    Complex q0=CM->ocenter;
    Complex q1=CM->d1;
    Complex q2=CM->d0;
    Complex q3=CM->od0;
    return HalfPlaneToDisk((PhiLast(Composition(q,q->Count,PhiFirst(z,q2,q1)),q3)),q0);
}

//---------------------------------------------------------------------------
Complex __fastcall GetInverseConformalMap(InverseConformalMapRec *ICM,Complex z)
{
    TList *q=ICM->data;
    Complex q0=ICM->ocenter;
    Complex q1=ICM->d1;
    Complex q2=ICM->d0;
    Complex q3=ICM->od0;
    Complex t1=DiskToHalfPlane(z,q0);
    Complex t2=PhiLastInv(t1,q3);
    Complex t3=InvComposition(q,q->Count,t2);
    Complex t4=PhiFirstInv(t3,q2,q1);
    return t4;
}

//---------------------------------------------------------------------------

TList *GenerateGrid(TList *lst,int polar,Complex &centr)
{
   TList *prb=new TList();
   TList *Zapolarray=Zapol(lst);
   double pntcount=Razmer(Zapolarray);
//   double pntcount=MaxValuePoint(lst);
   double xstep=2*pntcount/100;
   double ystep=2*pntcount/100;

   double hx=2*pntcount/400;
   double hy=2*pntcount/400;

   if(!polar)
   {
     for(int i=0;i<100;i++)
     {
        double x=xstep*i-pntcount;
        for(int j=0;j<400;j++)
        {
           double y=pntcount-hy*j;
           if(inside(x,y,Zapolarray))
           {
             Complex *t1=new Complex();
             t1->re=x;
             t1->im=y;
             prb->Add(t1);
           }
        }
     }
     for(int i=0;i<100;i++)
     {
        double y=ystep*i-pntcount;
        for(int j=0;j<400;j++)
        {
           double x=pntcount-hx*j;
           if(inside(x,y,Zapolarray))
           {
             Complex *t1=new Complex();
             t1->re=x;
             t1->im=y;
             prb->Add(t1);
           }
        }
     }
   }/* else
   {
       double pi2=2*3.14159265358979323846;
       double x,y;
       int l;
       double hx1=hx;
       for(int z=0;z<100;z++)
       {
          hx=(z+1)*hx1;
       for(int i=0;i<100;i++)
       {
           l=1;
           x=centr.re+hx*l*cos(i*pi2/1000);
           y=centr.im+hx*l*sin(i*pi2/1000);
           while(inside(x,y,lst))
           {
             l++;
             x=centr.re+hx*l*cos(i*pi2/1000);
             y=centr.im+hx*l*sin(i*pi2/1000);
           }
           l--;
           x=centr.re+hx*l*cos(i*pi2/1000);
           y=centr.im+hx*l*sin(i*pi2/1000);
           Complex *t1=new Complex();
           t1->re=x;
           t1->im=y;
           prb->Add(t1);

       }
       }

   }  */

   return prb;
}



//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormCreate(TObject *Sender)
{
   test=new TList();
}
//---------------------------------------------------------------------------
void __fastcall GeneratePoint(TList *lst,double x,double y)
{
   Complex *t1=new Complex();
   if(sign(x)==-1) t1->re=x;//+0.00001;
   else
   {
      if(x<1e-34) t1->re=x;//+0.00001;
      else t1->re=x;//-0.00001;
   }
   if(sign(y)==-1) t1->im=y;//+0.00001;
   else
   {
      if(y<1e-34) t1->im=y;//+0.00001;
      else t1->im=y;//-0.00001;
   }
   lst->Add(t1);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button1Click(TObject *Sender)
{
   for(int i=0;i<Form1->test->Count;i++) delete Form1->test->Items[i];
   Form1->test->Clear();

   //Êîò
   GeneratePoint(Form1->test,3,-4);
//   GeneratePoint(Form1->test,7,0);
//   GeneratePoint(Form1->test,7,5);
//   GeneratePoint(Form1->test,6,7);
//   GeneratePoint(Form1->test,6,8);
//   GeneratePoint(Form1->test,8,10);
//   GeneratePoint(Form1->test,7,11);
//   GeneratePoint(Form1->test,6,11);
//   GeneratePoint(Form1->test,8,11);
//   GeneratePoint(Form1->test,7,11);
//   GeneratePoint(Form1->test,5,9);
//   GeneratePoint(Form1->test,5,7);
//   GeneratePoint(Form1->test,6,6);
//   GeneratePoint(Form1->test,6,0);
   GeneratePoint(Form1->test,4,0);
   GeneratePoint(Form1->test,5,1);
   GeneratePoint(Form1->test,5,2);
   GeneratePoint(Form1->test,4,3);
   GeneratePoint(Form1->test,3,2);
   GeneratePoint(Form1->test,2,0);
   GeneratePoint(Form1->test,1,5);
   GeneratePoint(Form1->test,-1,7);
   GeneratePoint(Form1->test,-1,8);
   GeneratePoint(Form1->test,0,9);
   GeneratePoint(Form1->test,0,11);
   GeneratePoint(Form1->test,-1,12);
   GeneratePoint(Form1->test,-2,12);
   GeneratePoint(Form1->test,-3,11);
   GeneratePoint(Form1->test,-3,9);
   GeneratePoint(Form1->test,-2,8);
   GeneratePoint(Form1->test,-2,7);
   GeneratePoint(Form1->test,-6,3);
//   GeneratePoint(Form1->test,-6,-3);
   GeneratePoint(Form1->test,-6,-4);
//   GeneratePoint(Form1->test,-7,-3);
//   GeneratePoint(Form1->test,-7,-4);
   GeneratePoint(Form1->test,-5,-4);
   GeneratePoint(Form1->test,-5,2);
   GeneratePoint(Form1->test,-4,3);
   GeneratePoint(Form1->test,-2,1);
   GeneratePoint(Form1->test,-2,-1);
   GeneratePoint(Form1->test,-3,0);
   GeneratePoint(Form1->test,-4,-1);
//   GeneratePoint(Form1->test,-4,-2);
//   GeneratePoint(Form1->test,-3,-3);
//   GeneratePoint(Form1->test,-4,-3);
   GeneratePoint(Form1->test,-4,-4);

   Complex cntr;
   center(Form1->test,cntr),
   Form1->conmap=ConformalMap(Form1->test,cntr);
   Form1->conmap->grid=GenerateGrid(Form1->test,0,cntr);
   Form1->invconmap=GenerateInverseConformalMap(Form1->conmap->polygon,Form1->conmap->data,Form1->conmap->d0,Form1->conmap->d1,Form1->conmap->od0,Form1->conmap->ocenter,Form1->conmap->center);


   for(int i=0;i<Form1->test->Count;i++) delete Form1->test->Items[i];
   Form1->test->Clear();

   //Äèíîçàâð
   GeneratePoint(Form1->test,3,-1);
   GeneratePoint(Form1->test,1,2);
   GeneratePoint(Form1->test,1,4);
   GeneratePoint(Form1->test,3,2);
   GeneratePoint(Form1->test,4,2);
   GeneratePoint(Form1->test,2,6);
   GeneratePoint(Form1->test,-1,6);
   GeneratePoint(Form1->test,-2,4);
   GeneratePoint(Form1->test,-2,2);
   GeneratePoint(Form1->test,-3,3);
   GeneratePoint(Form1->test,-3,0);
   GeneratePoint(Form1->test,-2,1);
   GeneratePoint(Form1->test,-2,-3);
   GeneratePoint(Form1->test,-4,-2);
   GeneratePoint(Form1->test,-2,-4);
   GeneratePoint(Form1->test,1,-4);


   center(Form1->test,cntr),
   Form1->conmap1=ConformalMap(Form1->test,cntr);
   Form1->conmap1->grid=GenerateGrid(Form1->test,0,cntr);
   Form1->invconmap1=GenerateInverseConformalMap(Form1->conmap1->polygon,Form1->conmap1->data,Form1->conmap1->d0,Form1->conmap1->d1,Form1->conmap1->od0,Form1->conmap1->ocenter,Form1->conmap1->center);

   Button1->Enabled=false;
   Button2->Enabled=true;
   Button3->Enabled=true;
   Button4->Enabled=true;
   Button5->Enabled=true;
   Button6->Enabled=true;
   Button7->Enabled=true;
   ShowMessage("Ok");

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=4;
    double pntcount1=MaxValuePoint(conmap->polygon);
    double pntcount2=MaxValuePoint(conmap->polygon);
    int mash;
    if(pntcount1>pntcount2) mash=int(Image1->Height/pntcount1/2);
    else mash=int(Image1->Height/pntcount2/2);
    int offsetX=Image1->Width/2;
    int offsetY=Image1->Height/2;
    Image1->Canvas->MoveTo(offsetX,0);
    Image1->Canvas->LineTo(offsetX,Image1->Height);
    Image1->Canvas->MoveTo(0,offsetY);
    Image1->Canvas->LineTo(Image1->Width,offsetY);
    int l;
    l=conmap1->grid->Count;
    for(int i=0;i<l;i++)
    {
         Complex *t2=(Complex *)(conmap1->grid->Items[i]);
//===============================================================================================
//    Êîíôîðìíîå îòîáðàæåíèå âíóòð.äèíîçàâðà -> åä.êðóã -> âíóòð.êîòà
         Complex t1=GetInverseConformalMap(Form1->invconmap,GetConformalMap(Form1->conmap1,*t2));
//===============================================================================================
         int x=int(t1.re*mash+double(offsetX));
         int y=int(double(offsetY)-t1.im*mash);
         Image1->Canvas->Ellipse(x,y,x+2,y+2);
         Application->ProcessMessages();
    }

    l=conmap->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap->polygon->Items[i]);
       int x=int(t2->re*mash+double(offsetX));
       int y=int(double(offsetY)-t2->im*mash);
       Image1->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    Button2->Enabled=false;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button3Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=1;
    double pntcount1=MaxValuePoint(conmap->polygon);
    int mash=int(Image3->Height/pntcount1/2);
    int offsetX=Image3->Width/2;
    int offsetY=Image3->Height/2;
    Image3->Canvas->MoveTo(offsetX,0);
    Image3->Canvas->LineTo(offsetX,Image3->Height);
    Image3->Canvas->MoveTo(0,offsetY);
    Image3->Canvas->LineTo(Image3->Width,offsetY);
    int l=conmap->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap->polygon->Items[i]);
       Complex t1=*t2;
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image3->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    l=conmap->polygon->Count;
    for(int i=0;i<conmap->grid->Count;i++)
    {
       Complex *t2=(Complex *)(conmap->grid->Items[i]);
       Complex t1=*t2;
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image3->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    Button3->Enabled=false;
//    Button5->Enabled=true;

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button4Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=0;
    double pntcount1=MaxValuePoint(conmap1->polygon);
    int mash=int(Image2->Height/pntcount1/2);
    int offsetX=Image2->Width/2;
    int offsetY=Image2->Height/2;
    Image2->Canvas->MoveTo(offsetX,0);
    Image2->Canvas->LineTo(offsetX,Image2->Height);
    Image2->Canvas->MoveTo(0,offsetY);
    Image2->Canvas->LineTo(Image2->Width,offsetY);
    int l=conmap1->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap1->polygon->Items[i]);
       Complex t1=*t2;
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image2->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    l=conmap1->grid->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap1->grid->Items[i]);
       Complex t1=*t2;
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image2->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    Button4->Enabled=false;
//    Button3->Enabled=true;

}
//---------------------------------------------------------------------------



void __fastcall TForm1::Button5Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=2;
    double pntcount1=1;
    int mash=int(Image4->Height/pntcount1/2);
    int offsetX=Image4->Width/2;
    int offsetY=Image4->Height/2;
    Image4->Canvas->MoveTo(offsetX,0);
    Image4->Canvas->LineTo(offsetX,Image4->Height);
    Image4->Canvas->MoveTo(0,offsetY);
    Image4->Canvas->LineTo(Image4->Width,offsetY);
    int l=1000;
    for(int i=0;i<l;i++)
    {
       Complex t1=Complex(cos(2*3.14159265358979323846*i/1000),sin(2*3.14159265358979323846*i/1000));
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image4->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    l=conmap1->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap1->polygon->Items[i]);
       Complex t1=GetConformalMap(Form1->conmap1,*t2);
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image4->Canvas->Ellipse(x,y,x+4,y+4);
       Application->ProcessMessages();
    }
    l=conmap1->grid->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap1->grid->Items[i]);
       Complex t1=GetConformalMap(Form1->conmap1,*t2);
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image4->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    Button5->Enabled=false;
//    Button6->Enabled=true;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button6Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=3;
    double pntcount1=1;
    int mash=int(Image5->Height/pntcount1/2);
    int offsetX=Image5->Width/2;
    int offsetY=Image5->Height/2;
    Image5->Canvas->MoveTo(offsetX,0);
    Image5->Canvas->LineTo(offsetX,Image5->Height);
    Image5->Canvas->MoveTo(0,offsetY);
    Image5->Canvas->LineTo(Image5->Width,offsetY);
    int l=1000;
    for(int i=0;i<l;i++)
    {
       Complex t1=Complex(cos(2*3.14159265358979323846*i/1000),sin(2*3.14159265358979323846*i/1000));
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image5->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    l=conmap->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap->polygon->Items[i]);
       Complex t1=GetConformalMap(Form1->conmap,*t2);
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image5->Canvas->Ellipse(x,y,x+4,y+4);
       Application->ProcessMessages();
    }
    l=conmap->grid->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap->grid->Items[i]);
       Complex t1=GetConformalMap(Form1->conmap,*t2);
       int x=t1.re*mash+offsetX;
       int y=offsetY-t1.im*mash;
       Image5->Canvas->Ellipse(x,y,x+2,y+2);
       Application->ProcessMessages();
    }
    Button6->Enabled=false;
//    Button2->Enabled=true;
//    Button7->Enabled=true;

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button7Click(TObject *Sender)
{
    PageControl1->ActivePageIndex=5;
    double pntcount1=MaxValuePoint(conmap->polygon);
    double pntcount2=MaxValuePoint(conmap->polygon);
    int mash;
    if(pntcount1>pntcount2) mash=int(Image6->Height/pntcount1/2);
    else mash=int(Image6->Height/pntcount2/2);
    int offsetX=Image6->Width/2;
    int offsetY=Image6->Height/2;
    Image6->Canvas->MoveTo(offsetX,0);
    Image6->Canvas->LineTo(offsetX,Image6->Height);
    Image6->Canvas->MoveTo(0,offsetY);
    Image6->Canvas->LineTo(Image6->Width,offsetY);
    int l;
    l=conmap->grid->Count;
    for(int i=0;i<l;i++)
    {
         Complex *t2=(Complex *)(conmap->grid->Items[i]);
         Complex t1=GetInverseConformalMap(Form1->invconmap,GetConformalMap(Form1->conmap,*t2));
         int x=int(t1.re*mash+double(offsetX));
         int y=int(double(offsetY)-t1.im*mash);
         Image6->Canvas->Ellipse(x,y,x+2,y+2);
         Application->ProcessMessages();
    }

    l=conmap->polygon->Count;
    for(int i=0;i<l;i++)
    {
       Complex *t2=(Complex *)(conmap->polygon->Items[i]);
       for(int j=1999;j<2000;j++)
       {
         Complex t1=*t2;
         int x=int(t1.re*mash+double(offsetX));
         int y=int(double(offsetY)-t1.im*mash);
         Image6->Canvas->Ellipse(x,y,x+2,y+2);
         Application->ProcessMessages();
       }
    }
    Button7->Enabled=false;
}
//---------------------------------------------------------------------------

