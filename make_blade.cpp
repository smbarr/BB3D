#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <TColgp_Array1OfPnt.hxx>
#include <TColGeom_Array1OfBSplineCurve.hxx>
#include <gp_Pnt.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRep_Builder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <GeomLib_Interpolate.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <Geom_Surface.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <GeomFill_BoundWithSurf.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom2dAdaptor_HCurve.hxx>
#include <GeomAdaptor_HSurface.hxx>
#include <GeomAdaptor_HCurve.hxx>
#include <Adaptor3d_CurveOnSurface.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <OCC.h>

using namespace std;

double get_param_x(double, Handle(Geom_Curve));
gp_Pnt getPoint(std::string, double);
gp_Pnt rotate(gp_Pnt, double, double, double, double);
gp_Pnt rotate(gp_Pnt, double, double, double, double, bool);

int main(int argc, char **argv) {

  ifstream in("data.dat");
  std::string line;
  int station_number, nlofts=0, loft_num, nsections=0, nspars;
  double *aero_cent, *aero_orig, *twist, *chord, *z_station, dr, *spars;
  string *filename;

  IGESControl_Controller::Init();
  IGESControl_Writer iges_write;

  getline(in,line);
  stringstream ss;
  ss << line;
  ss >> nsections;
  getline(in,line);
  ss.clear();
  ss.str(std::string());
  ss.clear();
  ss << line;
  ss >> nlofts;
  getline(in,line);
  ss.clear();
  ss.str(std::string());
  ss.clear();
  ss << line;
  ss >> nspars;
 
  TColGeom_Array1OfBSplineCurve top_curves(0, nsections-1);
  TColGeom_Array1OfBSplineCurve bot_curves(0, nsections-1);
  aero_cent = new double[nsections];
  aero_orig = new double[nsections];
  twist = new double[nsections];
  chord = new double[nsections];
  z_station = new double[nsections];
  filename = new string[nsections];

  BRepOffsetAPI_ThruSections **top_lofts, **bot_lofts;
  top_lofts = new BRepOffsetAPI_ThruSections*[nlofts];
  bot_lofts = new BRepOffsetAPI_ThruSections*[nlofts];
  for(int n=0;n<nlofts;n++) {
    top_lofts[n] = new BRepOffsetAPI_ThruSections(Standard_False, Standard_False, 1.0e-2);
    bot_lofts[n] = new BRepOffsetAPI_ThruSections(Standard_False, Standard_False, 1.0e-2);
    top_lofts[n]->SetSmoothing(true);
    bot_lofts[n]->SetSmoothing(true);
  }
  BRepOffsetAPI_ThruSections **spar_surfs;
  spar_surfs = new BRepOffsetAPI_ThruSections*[nspars];
  for(int n=0;n<nspars;n++) {
    spar_surfs[n] = new BRepOffsetAPI_ThruSections(Standard_False, Standard_False, 1.0e-2);
    spar_surfs[n]->SetSmoothing(true);
  }
  BRep_Builder blade_builder;
  TopoDS_Compound blade_comp;
  blade_builder.MakeCompound(blade_comp);

  for(int n=0;n<nsections;n++) {
    getline(in,line);
    ss.clear();
    ss.str(std::string());
    ss << line;
    ss >> station_number;
    ss >> filename[station_number-1];
  }
  getline(in,line);

  for(int n=0;n<nsections;n++) {
    getline(in,line);
    ss.clear();
    ss.str(std::string());
    ss << line;
    ss >> station_number;
    ss >> z_station[station_number-1];
    ss >> twist[station_number-1];
    ss >> dr;
    ss >> chord[station_number-1];
    ss >> aero_cent[station_number-1];
    ss >> aero_orig[station_number-1];
    ss >> loft_num;
    //printf("%d %f %f %f %f %f %d\n",station_number,z_station[n],twist[n], chord[n],aero_cent[n],aero_orig[n],loft_num);
    cout << filename[station_number-1].c_str() << endl;
    ifstream airfoil(filename[station_number-1].c_str());
    int np = 0;
    while(!airfoil.eof()) {
      getline(airfoil,line);
      np++;
    }
    np = np - 1;
    airfoil.clear();
    airfoil.seekg(0, ios::beg);
    gp_Pnt pnt(0.0, 0.0, 0.0);
    gp_Pnt new_pnt(0.0, 0.0, 0.0);
    int loc, cnt1=0,cnt2=0;
    bool top_side = true, fore = false;
    char c;
    TColgp_Array1OfPnt top(0,(np-1)/2);
    TColgp_Array1OfPnt bot(0,(np-1)/2);
    for(int i=0;i<(np-1)/2;i++) {
      getline(airfoil,line);
      pnt = getPoint(line, z_station[station_number-1]);
      new_pnt = rotate(pnt, chord[station_number-1], twist[station_number-1], aero_cent[station_number-1], aero_orig[station_number-1]);
      top.SetValue(i, new_pnt);
    }
    getline(airfoil, line);
    pnt = getPoint(line, z_station[station_number-1]);
    new_pnt = rotate(pnt, chord[station_number-1], twist[station_number-1], aero_cent[station_number-1], aero_orig[station_number-1]);
    top.SetValue((np-1)/2, new_pnt);
    bot.SetValue(0, new_pnt);

    for(int i=1;i<(np-1)/2;i++) {
      getline(airfoil,line);
      pnt = getPoint(line, z_station[station_number-1]);
      new_pnt = rotate(pnt, chord[station_number-1], twist[station_number-1], aero_cent[station_number-1], aero_orig[station_number-1]);
      bot.SetValue(i, new_pnt);
    }
    getline(airfoil, line);
    pnt = getPoint(line, z_station[station_number-1]);
    new_pnt = rotate(pnt, chord[station_number-1], twist[station_number-1], aero_cent[station_number-1], aero_orig[station_number-1]);
    bot.SetValue((np-1)/2, new_pnt);

    airfoil.close();
    Standard_Boolean per = false;
    Standard_Real tol = 1.0e-4;
    GeomAPI_PointsToBSpline interp_top(top, 3, 3);
    GeomAPI_PointsToBSpline interp_bot(bot, 3, 3);
    Handle(Geom_BSplineCurve) top_curve = interp_top.Curve();
    Handle(Geom_BSplineCurve) bot_curve = interp_bot.Curve();
    top_curves.SetValue(station_number-1, top_curve);
    bot_curves.SetValue(station_number-1, bot_curve);

    BRepBuilderAPI_MakeEdge *te = new BRepBuilderAPI_MakeEdge(top_curve);
    BRepBuilderAPI_MakeEdge *be = new BRepBuilderAPI_MakeEdge(bot_curve);
    BRepBuilderAPI_MakeWire *tw = new BRepBuilderAPI_MakeWire(te->Edge());
    BRepBuilderAPI_MakeWire *bw = new BRepBuilderAPI_MakeWire(be->Edge());
    top_lofts[loft_num-1]->AddWire(tw->Wire());
    bot_lofts[loft_num-1]->AddWire(bw->Wire());

    iges_write.AddGeom(top_curve);
    iges_write.AddGeom(bot_curve);
  }

  getline(in, line);
  spars = new double[nspars];
  for(int n=0;n<nspars;n++) {
    getline(in, line);
    ss.clear();
    ss.str(std::string());
    ss << line;
    ss >> spars[n];
  }

  for(int n=0;n<nlofts;n++) {
    //iges_write.AddShape(top_lofts[n]->Shape());
    //iges_write.AddShape(bot_lofts[n]->Shape());
    blade_builder.Add(blade_comp,top_lofts[n]->Shape());
    blade_builder.Add(blade_comp,bot_lofts[n]->Shape());
    if(n<nlofts-1) {
      bspline_curve *le = new bspline_curve(3, 4);
      bspline_curve *tet = new bspline_curve(3, 4);
      bspline_curve *teb = new bspline_curve(3, 4);
      TopoDS_Face face1, face2, face3, face4;
      for(TopExp_Explorer Exe(top_lofts[n]->Shape(), TopAbs_FACE); Exe.More(); Exe.Next()) {
        face1 = TopoDS::Face(Exe.Current());
      }
      for(TopExp_Explorer Exe(top_lofts[n+1]->Shape(), TopAbs_FACE); Exe.More(); Exe.Next()) {
        face2 = TopoDS::Face(Exe.Current());
      }
      for(TopExp_Explorer Exe(bot_lofts[n+1]->Shape(), TopAbs_FACE); Exe.More(); Exe.Next()) {
        face3 = TopoDS::Face(Exe.Current());
      }
      for(TopExp_Explorer Exe(bot_lofts[n]->Shape(), TopAbs_FACE); Exe.More(); Exe.Next()) {
        face4 = TopoDS::Face(Exe.Current());
      }
      Handle(Geom_Surface) surf1 = BRep_Tool::Surface(face1);
      Handle(Geom_Surface) surf2 = BRep_Tool::Surface(face2);
      Handle(Geom_Surface) surf3 = BRep_Tool::Surface(face3);
      Handle(Geom_Surface) surf4 = BRep_Tool::Surface(face4);
      gp_Pnt sp, ep;
      gp_Vec svu, svv, evu, evv;
      Standard_Real u1,u2,v1,v2;

      surf1->Bounds(u1,u2,v1,v2);
      surf1->D1(u2,v2,sp,svu,svv);
      surf2->Bounds(u1,u2,v1,v2);
      surf2->D1(u2,v1,ep,evu,evv);
      double dz = (ep.Z()-sp.Z())/3.0;
      le->SetPoint(1,sp);
      le->SetPoint(2,sp.X()+dz*svv.X()/svv.Z(),sp.Y()+dz*svv.Y()/svv.Z(),sp.Z()+dz);
      le->SetPoint(3,ep.X()-dz*evv.X()/evv.Z(),ep.Y()-dz*evv.Y()/evv.Z(),ep.Z()-dz);
      le->SetPoint(4,ep);
      //iges_write.AddGeom(le->curve());

      surf1->Bounds(u1,u2,v1,v2);
      surf1->D1(u1,v2,sp,svu,svv);
      surf2->Bounds(u1,u2,v1,v2);
      surf2->D1(u1,v1,ep,evu,evv);
      dz = (ep.Z()-sp.Z())/3.0;
      tet->SetPoint(1,sp);
      tet->SetPoint(2,sp.X()+dz*svv.X()/svv.Z(),sp.Y()+dz*svv.Y()/svv.Z(),sp.Z()+dz);
      tet->SetPoint(3,ep.X()-dz*evv.X()/evv.Z(),ep.Y()-dz*evv.Y()/evv.Z(),ep.Z()-dz);
      tet->SetPoint(4,ep);
      //iges_write.AddGeom(tet->curve());

      surf4->Bounds(u1,u2,v1,v2);
      surf4->D1(u2,v2,sp,svu,svv);
      surf3->Bounds(u1,u2,v1,v2);
      surf3->D1(u2,v1,ep,evu,evv);
      dz = (ep.Z()-sp.Z())/3.0;
      teb->SetPoint(1,sp);
      teb->SetPoint(2,sp.X()+dz*svv.X()/svv.Z(),sp.Y()+dz*svv.Y()/svv.Z(),sp.Z()+dz);
      teb->SetPoint(3,ep.X()-dz*evv.X()/evv.Z(),ep.Y()-dz*evv.Y()/evv.Z(),ep.Z()-dz);
      teb->SetPoint(4,ep);
      //iges_write.AddGeom(teb->curve());

      surf1->Bounds(u1,u2,v1,v2);
      surf1->D1(u1,v2,sp,svu,svv);
      Handle(Geom_Curve) sct = surf1->VIso(v2);

      surf2->Bounds(u1,u2,v1,v2);
      surf2->D1(u1,v2,sp,svu,svv);
      Handle(Geom_Curve) ect = surf2->VIso(v1);

      surf3->Bounds(u1,u2,v1,v2);
      surf3->D1(u1,v2,sp,svu,svv);
      Handle(Geom_Curve) ecb = surf3->VIso(v1);

      surf4->Bounds(u1,u2,v1,v2);
      surf4->D1(u1,v2,sp,svu,svv);
      Handle(Geom_Curve) scb = surf4->VIso(v2);

      GeomAdaptor_Curve crv1(le->curve());
      GeomAdaptor_Curve crv2(tet->curve());
      GeomAdaptor_Curve crv3(teb->curve());
      Handle(GeomAdaptor_HCurve) hcurve1 = new GeomAdaptor_HCurve(crv1);
      Handle(GeomAdaptor_HCurve) hcurve2 = new GeomAdaptor_HCurve(crv2);
      Handle(GeomAdaptor_HCurve) hcurve3 = new GeomAdaptor_HCurve(crv3);
      
      GeomFill_SimpleBound *bnd1 = new GeomFill_SimpleBound(hcurve1,1e-2,0);
      GeomFill_SimpleBound *bnd2 = new GeomFill_SimpleBound(hcurve2,1e-2,0);
      GeomFill_SimpleBound *bnd3 = new GeomFill_SimpleBound(hcurve3,1e-2,0);
      GeomFill_ConstrainedFilling *top_surf = new GeomFill_ConstrainedFilling(3,250);
      GeomFill_ConstrainedFilling *bot_surf = new GeomFill_ConstrainedFilling(3,250);

      TopExp_Explorer exp1(face1, TopAbs_EDGE);
      exp1.Next(); exp1.Next();
      TopoDS_Edge ed1 = TopoDS::Edge(exp1.Current());

      TopExp_Explorer exp2(face2, TopAbs_EDGE);
      TopoDS_Edge ed2 = TopoDS::Edge(exp2.Current());

      TopExp_Explorer exp3(face3, TopAbs_EDGE);
      TopoDS_Edge ed3 = TopoDS::Edge(exp3.Current());

      TopExp_Explorer exp4(face4, TopAbs_EDGE);
      exp4.Next(); exp4.Next();
      TopoDS_Edge ed4 = TopoDS::Edge(exp4.Current());

      Standard_Real Tol3=1.0e-6, Tola=1.0e-3;
      surf1->Bounds(u1,u2,v1,v2);
      Handle(Geom2d_Curve) hcurve = BRep_Tool::CurveOnSurface(ed1,face1,v1,v2);
      hcurve = new Geom2d_TrimmedCurve(hcurve,v1,v2);
      Handle(Geom2d_Curve) myParamCurve = hcurve;
      Handle(GeomAdaptor_HSurface) Surface = new GeomAdaptor_HSurface(surf1);
      Handle(Geom2dAdaptor_HCurve) ParamCurve = new Geom2dAdaptor_HCurve(myParamCurve);
      Adaptor3d_CurveOnSurface CurveOnSurf = Adaptor3d_CurveOnSurface(ParamCurve,Surface);
      GeomFill_BoundWithSurf *myBoundary1 = new GeomFill_BoundWithSurf(CurveOnSurf, Tol3, Tola);

      surf2->Bounds(u1,u2,v1,v2);
      hcurve = BRep_Tool::CurveOnSurface(ed2,face2,v1,v2);
      hcurve = new Geom2d_TrimmedCurve(hcurve,v1,v2);
      myParamCurve = hcurve;
      Surface = new GeomAdaptor_HSurface(surf2);
      ParamCurve = new Geom2dAdaptor_HCurve(myParamCurve);
      CurveOnSurf = Adaptor3d_CurveOnSurface(ParamCurve,Surface);
      GeomFill_BoundWithSurf *myBoundary2 = new GeomFill_BoundWithSurf(CurveOnSurf, Tol3, Tola);

      surf3->Bounds(u1,u2,v1,v2);
      hcurve = BRep_Tool::CurveOnSurface(ed3,face3,v1,v2);
      hcurve = new Geom2d_TrimmedCurve(hcurve,v1,v2);
      myParamCurve = hcurve;
      Surface = new GeomAdaptor_HSurface(surf3);
      ParamCurve = new Geom2dAdaptor_HCurve(myParamCurve);
      CurveOnSurf = Adaptor3d_CurveOnSurface(ParamCurve,Surface);
      GeomFill_BoundWithSurf *myBoundary3 = new GeomFill_BoundWithSurf(CurveOnSurf, Tol3, Tola);

      surf4->Bounds(u1,u2,v1,v2);
      hcurve = BRep_Tool::CurveOnSurface(ed4,face4,v1,v2);
      hcurve = new Geom2d_TrimmedCurve(hcurve,v1,v2);
      myParamCurve = hcurve;
      Surface = new GeomAdaptor_HSurface(surf4);
      ParamCurve = new Geom2dAdaptor_HCurve(myParamCurve);
      CurveOnSurf = Adaptor3d_CurveOnSurface(ParamCurve,Surface);
      GeomFill_BoundWithSurf *myBoundary4 = new GeomFill_BoundWithSurf(CurveOnSurf, Tol3, Tola);
      
      top_surf->Init(myBoundary1,bnd1,myBoundary2,bnd2,Standard_True);
      BRepBuilderAPI_MakeFace *top_face = new BRepBuilderAPI_MakeFace(top_surf->Surface(),1.0e-6);
      //iges_write.AddShape(top_face->Face());

      bot_surf->Init(myBoundary3,bnd1,myBoundary4,bnd3,Standard_True);
      BRepBuilderAPI_MakeFace *bot_face = new BRepBuilderAPI_MakeFace(bot_surf->Surface(),1.0e-6);
      //iges_write.AddShape(bot_face->Face());

      blade_builder.Add(blade_comp, top_face->Shape());
      blade_builder.Add(blade_comp, bot_face->Shape());
    }
  }

  for(int n=0;n<nspars;n++) {
    gp_Pnt tp(spars[n], 10.0, z_station[0]-2.0);
    gp_Pnt bp(spars[n], -10.0, z_station[0]-2.0);
    gp_Pnt tpr, bpr;
    tpr = rotate(tp, chord[0], twist[0], aero_cent[0], aero_orig[0], false);
    bpr = rotate(bp, chord[0], twist[0], aero_cent[0], aero_orig[0], false);
    BRepBuilderAPI_MakeEdge *e = new BRepBuilderAPI_MakeEdge(tpr, bpr);
    BRepBuilderAPI_MakeWire *w = new BRepBuilderAPI_MakeWire(e->Edge());
    spar_surfs[n]->AddWire(w->Wire());

    for(int s=0;s<nsections;s++) {
      tp = gp_Pnt(spars[n], 10.0, z_station[s]);
      bp = gp_Pnt(spars[n], -10.0, z_station[s]);
      tpr = rotate(tp, chord[s], twist[s], aero_cent[s], aero_orig[s], false);
      bpr = rotate(bp, chord[s], twist[s], aero_cent[s], aero_orig[s], false);
      e = new BRepBuilderAPI_MakeEdge(tpr, bpr);
      w = new BRepBuilderAPI_MakeWire(e->Edge());
      spar_surfs[n]->AddWire(w->Wire());
    }

    tp = gp_Pnt(spars[n], 10.0, z_station[nsections-1]+2.0);
    bp = gp_Pnt(spars[n], -10.0, z_station[nsections-1]+2.0);
    tpr = rotate(tp, chord[nsections-1], twist[nsections-1], aero_cent[nsections-1], aero_orig[nsections-1], false);
    bpr = rotate(bp, chord[nsections-1], twist[nsections-1], aero_cent[nsections-1], aero_orig[nsections-1], false);
    e = new BRepBuilderAPI_MakeEdge(tpr, bpr);
    w = new BRepBuilderAPI_MakeWire(e->Edge());
    spar_surfs[n]->AddWire(w->Wire());
  }

  TopoDS_Shape blade = blade_comp;
  for(int n=0;n<nspars;n++) {
    BRepAlgoAPI_Cut *cut1 = new BRepAlgoAPI_Cut(blade, spar_surfs[n]->Shape());
    blade = cut1->Shape();
  }

  iges_write.AddShape(blade);

  iges_write.Write("blade.igs");
  
  return 0;
}

gp_Pnt getPoint(std::string line, double z) {
  gp_Pnt pnt;
  stringstream ss;

  double x,y;
  ss << line;
  ss >> x;
  ss >> y;
  pnt.SetX(x);
  pnt.SetY(y);
  pnt.SetZ(z);

  return pnt;
}

gp_Pnt rotate(gp_Pnt pnt, double factor, double twist, double aero_cent, double aero_orig) {
  double PI = 3.141592653;
  double tr = twist*PI/180.0;
  gp_Pnt new_pnt(0.0, 0.0, 0.0);

  new_pnt.SetX(pnt.X() - (aero_orig + (0.25-aero_cent)));
  new_pnt.SetX(new_pnt.X()*factor);
  new_pnt.SetY(pnt.Y()*factor);

  new_pnt.SetX((cos(tr)*new_pnt.X()) - (sin(tr)*new_pnt.Y()));
  new_pnt.SetY((sin(tr)*new_pnt.X()) + (cos(tr)*new_pnt.Y()));
  new_pnt.SetZ(pnt.Z());

  return new_pnt;
}

gp_Pnt rotate(gp_Pnt pnt, double factor, double twist, double aero_cent, double aero_orig, bool mult) {
  double PI = 3.141592653;
  double tr = twist*PI/180.0;
  gp_Pnt new_pnt(0.0, 0.0, 0.0);

  new_pnt.SetX(pnt.X() - (aero_orig + (0.25-aero_cent)));
  new_pnt.SetX(new_pnt.X()*factor);
  if(mult == true) {
    new_pnt.SetY(pnt.Y()*factor);
  } else {
    new_pnt.SetY(pnt.Y());
  }

  new_pnt.SetX((cos(tr)*new_pnt.X()) - (sin(tr)*new_pnt.Y()));
  new_pnt.SetY((sin(tr)*new_pnt.X()) + (cos(tr)*new_pnt.Y()));
  new_pnt.SetZ(pnt.Z());

  return new_pnt;
}

double get_param_x(double x, Handle(Geom_Curve) c) {
  double px, f1, f2, error = 1.0e6;
  gp_Pnt pnt;
  px = 0.0;
  do {
    c->D0(px,pnt);
    f1 = pnt.X()-x;
    c->D0(px+1.0e-3,pnt);
    f2 = pnt.X()-x;
    px = px - (f1/((f2-f1)/1.0e-3));
    error = abs(f1);
  } while(error > 1.0e-6);

  return px;
}
