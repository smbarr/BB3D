#ifndef OCC_H
#define OCC_H
//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL); // VTK was built with vtkRenderingOpenGL2
//VTK_MODULE_INIT(vtkInteractionStyle);
#include <math.h>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_HCurve.hxx>
#include <GeomFill_SimpleBound.hxx>
#include <GeomFill_ConstrainedFilling.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
//#include <IVtkOCC_Shape.hxx>
//#include <IVtkTools_ShapeDataSource.hxx>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <IVtkTools_DisplayModeFilter.hxx>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkInteractorStyleTrackballActor.h>
//#include <vtkObjectFactory.h>
//#include <vtkMatrix4x4.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkSmoothPolyDataFilter.h>
//#include <vtkPolyDataNormals.h>

class bezier_curve {
  int npnts;

  public:

    //int check() {
    //}

    bezier_curve(int n) {
      npnts = n;
      pnts = new double*[npnts];
      for(int i=0;i<npnts;i++) {
        pnts[i] = new double[3];
        pnts[i][0] = 0.0;
        pnts[i][1] = 0.0;
        pnts[i][2] = 0.0;
      }
    }

    void SetPoint(int n, double x, double y, double z) {
      pnts[n-1][0] = x;
      pnts[n-1][1] = y;
      pnts[n-1][2] = z;
    }

    void SetPoint(int n, gp_Pnt p) {
      pnts[n-1][0] = p.X();
      pnts[n-1][1] = p.Y();
      pnts[n-1][2] = p.Z();
    }

    Handle(Geom_BezierCurve) curve() {
      TColgp_Array1OfPnt tcps(0,npnts-1);
      for(int n=0;n<npnts;n++) {
        tcps.SetValue(n,gp_Pnt(pnts[n][0],pnts[n][1],pnts[n][2]));
      }
      Handle(Geom_BezierCurve) c = new Geom_BezierCurve(tcps);

      return c;
    };

  private:
    double **pnts;

};

class bspline_curve {
  int npnts, degree;

  public:

    //int check() {
    //}

    bspline_curve(int d, int n) {
      degree = d;
      npnts = n;
      wght = new double[npnts];
      pnts = new double*[npnts];
      for(int i=0;i<npnts;i++) {
        wght[i] = 1.0;
        pnts[i] = new double[3];
        pnts[i][0] = 0.0;
        pnts[i][1] = 0.0;
        pnts[i][2] = 0.0;
      }
    }

    void SetWeight(int n, double w) {
      wght[n-1] = w;
    }

    void SetPoint(int n, gp_Pnt pnt) {
      pnts[n-1][0] = pnt.X();
      pnts[n-1][1] = pnt.Y();
      pnts[n-1][2] = pnt.Z();
    }

    void SetPoint(int n, double x, double y, double z) {
      pnts[n-1][0] = x;
      pnts[n-1][1] = y;
      pnts[n-1][2] = z;
    }

    Handle(Geom_BSplineCurve) curve() {
      Handle(Geom_BSplineCurve) c;
      TColgp_Array1OfPnt tcps(0,npnts-1);
      TColStd_Array1OfReal w(0,npnts-1);

      for(int n=0;n<npnts;n++) {
        tcps.SetValue(n,gp_Pnt(pnts[n][0],pnts[n][1],pnts[n][2]));
        w.SetValue(n,wght[n]);
      }

      int knot_size = npnts+degree+1-(2*(degree));
      TColStd_Array1OfReal k(0,knot_size-1);
      TColStd_Array1OfInteger m(0,knot_size-1);
      for(int n=0;n<knot_size;n++) {
        k.SetValue(n,n+1.0);
        m.SetValue(n,1);
      }
      m.SetValue(0,degree+1);
      m.SetValue(knot_size-1,degree+1);
      c = new Geom_BSplineCurve(tcps, w, k, m, degree);

      return c;
    };

  private:
    double *wght;
    double **pnts;

};

//class nurbsSurface {
//  public:
//    surface(
//};

class ray {
  public:
    gp_Pnt point;
    gp_Dir direction;
    ray(gp_Pnt pnt, gp_Dir dir) {
      point.SetX(pnt.X());
      point.SetY(pnt.Y());
      point.SetZ(pnt.Z());
      direction.SetX(dir.X());
      direction.SetY(dir.Y());
      direction.SetZ(dir.Z());
    }
};

class bsplineSurface {

  public:
    Standard_Integer u_degree;
    Standard_Integer v_degree;
    Standard_Integer num_up;
    Standard_Integer num_vp;
    Standard_Integer num_uk;
    Standard_Integer num_vk;
    bool done;

    bsplineSurface() {
      done = false;
    }

    void constrainedFilling(Handle(Geom_Curve) c1, Handle(Geom_Curve) c2, Handle(Geom_Curve) c3, Handle(Geom_Curve) c4, int degree, int num_pnts) {
      GeomAdaptor_Curve crv1(c1);
      GeomAdaptor_Curve crv2(c2);
      GeomAdaptor_Curve crv3(c3);
      GeomAdaptor_Curve crv4(c4);
      Handle(GeomAdaptor_HCurve) hcurve1 = new GeomAdaptor_HCurve(crv1);
      Handle(GeomAdaptor_HCurve) hcurve2 = new GeomAdaptor_HCurve(crv2);
      Handle(GeomAdaptor_HCurve) hcurve3 = new GeomAdaptor_HCurve(crv3);
      Handle(GeomAdaptor_HCurve) hcurve4 = new GeomAdaptor_HCurve(crv4);

      GeomFill_SimpleBound *bnd1 = new GeomFill_SimpleBound(hcurve1,1e-6,1.0e-2);
      GeomFill_SimpleBound *bnd2 = new GeomFill_SimpleBound(hcurve2,1e-6,1.0e-2);
      GeomFill_SimpleBound *bnd3 = new GeomFill_SimpleBound(hcurve3,1e-6,1.0e-2);
      GeomFill_SimpleBound *bnd4 = new GeomFill_SimpleBound(hcurve4,1e-6,1.0e-2);
      GeomFill_ConstrainedFilling *surf = new GeomFill_ConstrainedFilling(degree, num_pnts);
      surf->Init(bnd1,bnd2,bnd3,bnd4,Standard_True);

      u_degree = surf->Surface()->UDegree();
      v_degree = surf->Surface()->VDegree();
      num_up = surf->Surface()->NbUPoles();
      num_vp = surf->Surface()->NbVPoles();
      num_uk = surf->Surface()->NbUKnots();
      num_vk = surf->Surface()->NbVKnots();
      wght = new double*[num_up];
      pnts = new double**[num_up];
      for(int i=1;i<=num_up;i++) {
        wght[i-1] = new double[num_vp];
        pnts[i-1] = new double*[num_vp];
        for(int j=1;j<=num_vp;j++) {
          gp_Pnt p = gp_Pnt(surf->Surface()->Pole(i,j));
          pnts[i-1][j-1] = new double[3];
          pnts[i-1][j-1][0] = p.X();
          pnts[i-1][j-1][1] = p.Y();
          pnts[i-1][j-1][2] = p.Z();
          wght[i-1][j-1] = surf->Surface()->Weight(i,j);
        }
      }
      
      done = true;
    }

    void extrude(Handle(Geom_BSplineCurve) c, gp_Dir d, double dist) {
      u_degree = c->Degree();
      num_up = c->NbPoles();
      num_uk = c->NbKnots();
      v_degree = 1;
      num_vp = 2;
      num_vk = 2;
      wght = new double*[num_up];
      pnts = new double**[num_up];
      for(int i=1;i<=num_up;i++) {
        wght[i-1] = new double[num_vp];
        pnts[i-1] = new double*[num_vp];
        for(int j=1;j<=num_vp;j++) {
          pnts[i-1][j-1] = new double[3];
          wght[i-1][j-1] = c->Weight(i);
        }
      }
      for(int i=1;i<=num_up;i++) {
        pnts[i-1][0][0] = c->Pole(i).X();
        pnts[i-1][0][1] = c->Pole(i).Y();
        pnts[i-1][0][2] = c->Pole(i).Z();
        pnts[i-1][1][0] = c->Pole(i).X()+(dist*d.X());
        pnts[i-1][1][1] = c->Pole(i).Y()+(dist*d.Y());
        pnts[i-1][1][2] = c->Pole(i).Z()+(dist*d.Z());
      }
      done = true;
    }

    void railRevolve(Handle(Geom_BSplineCurve) profile, Handle(Geom_BSplineCurve) rail, ray r) {
      // axis is a 1-d array with three elements a,b, and c, where a can be 0, 1, or 2 referring 
      //   to x, y, and z axes, respectively. b and c is the value of the axis. 
      //   Ex: axis = [0, 3.5, 4.5] would refer to an axis along x at y=3.5, z=4.5.

      // Requirement: Profile must pierce endpoint of rail at one end, and the axis with the other end

      //  1) Rotate profile pole point about axis
      //  2) Scale resultant point about vector in direction from origin to pole point on rail
      
      u_degree = rail->Degree();
      v_degree = profile->Degree();
      num_up = rail->NbPoles();
      num_vp = profile->NbPoles();
      num_uk = rail->NbKnots();
      num_vk = profile->NbKnots();

      int flip = 0;
      if(rail->Pole(num_up).IsEqual(profile->Pole(1), 1.0e-6)) {
        rail->Reverse();
        flip = 1;
      } else if(rail->Pole(1).IsEqual(profile->Pole(num_vp), 1.0e-6)) {
        profile->Reverse();
        flip = 2;
      } else if(rail->Pole(num_up).IsEqual(profile->Pole(num_vp), 1.0e-6)) {
        rail->Reverse();
        profile->Reverse();
        flip = 3;
      }

      pnts = new double**[num_up];
      wght = new double*[num_up];
      for(int i=0;i<num_up;i++) {
        pnts[i] = new double*[num_vp];
        wght[i] = new double[num_vp];
        for(int j=0;j<num_vp;j++) {
          wght[i][j] = 1.0;
          pnts[i][j] = new double[3];
        }
      }

      double theta = 0,umag;
      double *angles;
      gp_Pnt pnt_rot,pnt_scaled;
      gp_Vec u;

      angles = new double[num_up];

      if(r.direction.X() == 1.0) {
        for(int i=1;i<=num_up;i++) {
          double x = rail->Pole(1).X();
          gp_Vec vec1 = gp_Vec(rail->Pole(1),gp_Pnt(x,r.point.Y(),r.point.Z()));
          gp_Vec vec2 = gp_Vec(rail->Pole(i),gp_Pnt(x,r.point.Y(),r.point.Z()));
          double theta = vec2.Angle(vec1);
          angles[i-1] = theta*-1.0;
        }

        u.SetX(r.direction.X());
        u.SetY(r.direction.Y());
        u.SetZ(r.direction.Z());

        for(int i=1;i<=num_up;i++) {
          for(int j=1;j<=num_vp;j++) {
            pnt_rot = rotate(profile->Pole(j), u, angles[i-1], r.point);
            SetPoint(i-1, j-1, pnt_rot);
          }
        }

        gp_Pnt pnt1,pnt2,pnt3;
        for(int i=2;i<=num_up;i++) {
          pnt1 = gp_Pnt(pnts[i-1][0][0], r.point.Y(), r.point.Z());
          pnt2 = gp_Pnt(pnts[i-1][0][0],pnts[i-1][0][1],pnts[i-1][0][2]);
          pnt3 = rail->Pole(i);
          gp_Vec svec1 = gp_Vec(pnt1,pnt2);
          gp_Vec svec2 = gp_Vec(pnt1,pnt3);
          double scale_y, scale_z;
          double new_y, new_z;
          scale_y = svec2.Y()/svec1.Y();
          scale_z = svec2.Z()/svec1.Z();
          for(int j=1;j<=num_vp;j++) {
            new_y = ((pnts[i-1][j-1][1]-r.point.Y()) * scale_y) + r.point.Y();
            new_z = ((pnts[i-1][j-1][2]-r.point.Z()) * scale_z) + r.point.Z();
            SetPoint(i-1, j-1, gp_Pnt(pnts[i-1][j-1][0],new_y,new_z));
          }
        }
      }

      if(flip == 1) {
        rail->Reverse();
      } else if(flip == 2) {
        profile->Reverse();
      } else if(flip == 3) {
        rail->Reverse();
        profile->Reverse();
      }
      done = true;
    }

    void SetPoint(int i, int j, gp_Pnt p) {
      pnts[i][j][0] = p.X();
      pnts[i][j][1] = p.Y();
      pnts[i][j][2] = p.Z();
    }

    void PrintPoint(gp_Pnt p) {
      printf("%f %f %f\n",p.X(),p.Y(),p.Z());
    }

    void PrintPoint(int i, int j) {
      printf("%f %f %f\n",pnts[i][j][0],pnts[i][j][1],pnts[i][j][2]);
    }

    Handle(Geom_BSplineSurface) surface() {
      if(done == true) {
        int uknot_size = num_up+u_degree+1-(2*u_degree);
        int vknot_size = num_vp+v_degree+1-(2*v_degree);
        Handle(Geom_BSplineSurface) s;
        TColgp_Array2OfPnt tcps(0, num_up-1, 0, num_vp-1);
        TColStd_Array2OfReal w(0, num_up-1, 0, num_vp-1);
        TColStd_Array1OfReal uk(0, uknot_size-1);
        TColStd_Array1OfReal vk(0, vknot_size-1);
        TColStd_Array1OfInteger um(0, uknot_size-1);
        TColStd_Array1OfInteger vm(0, vknot_size-1);

        for(int i=0;i<num_up;i++) {
          for(int j=0;j<num_vp;j++) {
            tcps.SetValue(i,j,gp_Pnt(pnts[i][j][0],pnts[i][j][1],pnts[i][j][2]));
            w.SetValue(i,j,wght[i][j]);
          }
        }
        for(int i=0;i<uknot_size;i++) {
          uk.SetValue(i,i+1.0);
          um.SetValue(i,1);
        }
        for(int j=0;j<vknot_size;j++) {
          vk.SetValue(j,j+1.0);
          vm.SetValue(j,1);
        }
        um.SetValue(0, u_degree+1);
        um.SetValue(uknot_size-1, u_degree+1);
        vm.SetValue(0, v_degree+1);
        vm.SetValue(vknot_size-1, v_degree+1);

        s = new Geom_BSplineSurface(tcps, w, uk, vk, um, vm, u_degree, v_degree);

        return s;
      } else {
        cout << "Surface not built." << endl;
        exit(1);
      }
    }

    //TopoDS_Face face() {
    //  BRepBuilderAPI_MakeFace *mf = new BRepBuilderAPI_MakeFace(surface(), 1.0e-6);
    //  return mf->Shape();
    //}

    gp_Pnt rotate(gp_Pnt pnt, gp_Vec u, double theta, gp_Pnt trans) {
      double umag;
      gp_Pnt pt_trans,pt_new;
      double R[3][3];

      pt_trans.SetX(pnt.X()-trans.X());
      pt_trans.SetY(pnt.Y()-trans.Y());
      pt_trans.SetZ(pnt.Z()-trans.Z());

      umag = sqrt(pow(u.X(),2)+pow(u.Y(),2)+pow(u.Z(),2));
      u.SetX(u.X()/umag);
      u.SetY(u.Y()/umag);
      u.SetZ(u.Z()/umag);

      R[0][0] = cos(theta)+(pow(u.X(),2)*(1.0-cos(theta)));
      R[0][1] = (u.X()*u.Y()*(1.0-cos(theta)))-(u.Z()*sin(theta));
      R[0][2] = (u.X()*u.Z()*(1.0-cos(theta)))-(u.Y()*sin(theta));
      R[1][0] = (u.Y()*u.X()*(1.0-cos(theta)))-(u.Z()*sin(theta));
      R[1][1] = cos(theta)+(pow(u.Y(),2)*(1.0-cos(theta)));
      R[1][2] = (u.Y()*u.Z()*(1.0-cos(theta)))-(u.X()*sin(theta));
      R[2][0] = (u.Z()*u.X()*(1.0-cos(theta)))-(u.Y()*sin(theta));
      R[2][1] = (u.Z()*u.Y()*(1.0-cos(theta)))-(u.X()*sin(theta));
      R[2][2] = cos(theta)+(pow(u.Z(),2)*(1.0-cos(theta)));

      pt_new = matrix_mult(pt_trans,R);
      pt_new.SetX(pt_new.X()+trans.X());
      pt_new.SetY(pt_new.Y()+trans.Y());
      pt_new.SetZ(pt_new.Z()+trans.Z());
      return pt_new;
      
    }

    gp_Pnt matrix_mult(gp_Pnt a, double b[][3]) {
      double ca[3],aa[3];
      gp_Pnt c;

      aa[0] = a.X();
      aa[1] = a.Y();
      aa[2] = a.Z();

      for(int i=0;i<3;i++) {
        ca[i] = 0.0;
        for(int j=0;j<3;j++) {
          ca[i] = ca[i] + (b[i][j]*aa[j]);
        }
      }
      c.SetX(ca[0]);
      c.SetY(ca[1]);
      c.SetZ(ca[2]);

      return c;
    }

  private:
    double **wght;
    double ***pnts;
};

//class display {
//  public:
//    display() {
//      renderer = vtkSmartPointer<vtkRenderer>::New();
//      renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//      renderWindow->SetSize(1000, 1000);
//      renderWindow->AddRenderer(renderer);
//      renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//      renderWindowInteractor->SetRenderWindow(renderWindow);
//    }
//
//    void startDisplay() {
//      renderer->GradientBackgroundOn();
//      renderer->SetBackground(1,1,1);
//      renderer->SetBackground2(0,0,0);
//      renderWindow->Render();
//      vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
//      renderWindowInteractor->SetInteractorStyle( style );
//      renderWindowInteractor->Start();
//    }
//
//    void addShape(TopoDS_Shape shape) {
//      IVtkOCC_Shape::Handle aShapeImpl = new IVtkOCC_Shape(shape);
//      vtkSmartPointer<IVtkTools_ShapeDataSource> DS = vtkSmartPointer<IVtkTools_ShapeDataSource>::New();
//      DS->SetShape(aShapeImpl);
//
//      //vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//      //mapper->SetInputConnection(DS->GetOutputPort());
//
//      vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
//        vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
//      smoothFilter->SetInputConnection(DS->GetOutputPort());
//      smoothFilter->SetNumberOfIterations(15);
//      smoothFilter->SetRelaxationFactor(0.1);
//      smoothFilter->FeatureEdgeSmoothingOff();
//      smoothFilter->BoundarySmoothingOn();
//      smoothFilter->Update();
//
//      // Update normals on newly smoothed polydata
//      vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
//      normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
//      normalGenerator->ComputePointNormalsOn();
//      normalGenerator->ComputeCellNormalsOn();
//      normalGenerator->Update();
//
//      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//      mapper->SetInputConnection(normalGenerator->GetOutputPort());
//
//      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//      actor->SetMapper(mapper);
//      actor->GetProperty()->SetColor(0.67, 0.67, 0.67); //(R,G,B)
//      //actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); //(R,G,B)
//      //actor->GetProperty()->EdgeVisibilityOn();
//
//      renderer->AddActor(actor);
//    }
//
//    void addShape(TopoDS_Shape shape, double r, double g, double b) {
//      IVtkOCC_Shape::Handle aShapeImpl = new IVtkOCC_Shape(shape);
//      vtkSmartPointer<IVtkTools_ShapeDataSource> DS = vtkSmartPointer<IVtkTools_ShapeDataSource>::New();
//      DS->SetShape(aShapeImpl);
//
//      //vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//      //mapper->SetInputConnection(DS->GetOutputPort());
//
//      vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
//        vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
//      smoothFilter->SetInputConnection(DS->GetOutputPort());
//      smoothFilter->SetNumberOfIterations(15);
//      smoothFilter->SetRelaxationFactor(0.1);
//      smoothFilter->FeatureEdgeSmoothingOff();
//      smoothFilter->BoundarySmoothingOn();
//      smoothFilter->Update();
//
//      // Update normals on newly smoothed polydata
//      vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
//      normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
//      normalGenerator->ComputePointNormalsOn();
//      normalGenerator->ComputeCellNormalsOn();
//      normalGenerator->Update();
//
//      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//      mapper->SetInputConnection(normalGenerator->GetOutputPort());
//
//      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//      actor->SetMapper(mapper);
//      actor->GetProperty()->SetColor(r, g, b); //(R,G,B)
//      //actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); //(R,G,B)
//      //actor->GetProperty()->EdgeVisibilityOn();
//
//      renderer->AddActor(actor);
//    }
//  private:
//    vtkSmartPointer<vtkRenderer> renderer;
//    vtkSmartPointer<vtkRenderWindow> renderWindow;
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
//    //vtkSmartPointer<MyInteractorStyle> style;
//};

#endif
