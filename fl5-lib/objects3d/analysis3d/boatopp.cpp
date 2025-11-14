/****************************************************************************

    flow5 application
    Copyright (C) 2025 André Deperrois 
    
    This file is part of flow5.

    flow5 is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    flow5 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with flow5.
    If not, see <https://www.gnu.org/licenses/>.


*****************************************************************************/



#include <format>


#include <api/boatopp.h>


#include <api/boat.h>
#include <api/boatpolar.h>
#include <api/sail.h>
#include <api/sailglobals.h>
#include <api/units.h>
#include <api/utils.h>


BoatOpp::BoatOpp()
{
    m_BoatName.clear();
    m_BtPolarName.clear();
    m_nPanel3 = 0;
    m_nPanel4 = 0;

    m_AnalysisMethod = xfl::VLM2;

    m_bGround           = false;
    m_bIgnoreBodyPanels = false;
    m_bThinSurfaces     = false;
    m_bTrefftz          = false;

    m_bGround      = true;
    m_GroundHeight = 0.0;

    m_TWS_inf = m_TWA_inf = 0.0;

    m_NodeValMin = m_NodeValMax = 0.0;

    m_SailAngle.clear();
}


BoatOpp::BoatOpp(Boat *pBoat, BoatPolar *pBtPolar, int nPanel3, int nPanel4)
{
    m_BoatName = pBoat->name();
    m_BtPolarName = pBtPolar->name();
    m_nPanel3 = nPanel3;
    m_nPanel4 = nPanel4;

    m_AnalysisMethod = pBtPolar->analysisMethod();

    m_bGround = pBtPolar->m_bGround;
    m_bIgnoreBodyPanels = pBtPolar->bIgnoreBodyPanels();
//    m_bThinSurfaces = pBtPolar->m_bThinSurfaces;
    m_bTrefftz= pBtPolar->m_bTrefftz;

    m_GroundHeight      = 0.0;
    m_QInf        = 0.0;
    m_Beta        = 0.0;
    m_Phi         = 0.0;
    m_Ctrl        = 0.0;

    m_SailAngle.resize(pBtPolar->sailAngleSize());
    for(uint is=0;is<m_SailAngle.size(); is++) m_SailAngle[is]=0.0;

    m_SailForceFF.resize(pBoat->nSails());
    m_SailForceSum.resize(pBoat->nSails());
}


bool BoatOpp::serializeBoatOppFl5(QDataStream &ar, bool bIsStoring)
{
    // 100001: first file format
    // 100002: added lift and drag
    // 100003: added vortons and negating vortices
    // 100004: Modified the format of AeroForces serialization
    // 100005: beta20 - Added the roration about Ry

    int ArchiveFormat=100005;

    float f(0),g(0),h(0);

    int n(0);
    double dble(0);
    int nIntSpares(0);
    int nDbleSpares(0);
    QString strange;

    if(bIsStoring)
    {
        ar << ArchiveFormat;

        ar << QString::fromStdString(m_BoatName);
        ar << QString::fromStdString(m_BtPolarName);

        m_theStyle.serializeFl5(ar, bIsStoring);

        //ANALYSIS METHOD
        if     (m_AnalysisMethod==xfl::LLT)        ar<<1;
        else if(m_AnalysisMethod==xfl::VLM1)       ar<<2;
        else if(m_AnalysisMethod==xfl::VLM1)       ar<<3;
        else if(m_AnalysisMethod==xfl::QUADS)      ar<<4;
        else if(m_AnalysisMethod==xfl::TRILINEAR)  ar<<5;
        else if(m_AnalysisMethod==xfl::TRIUNIFORM) ar<<6;
        else                                       ar<<0;

        ar << m_bThinSurfaces;
        ar << m_bTrefftz;

        ar << m_bIgnoreBodyPanels;

        ar << m_GroundHeight << m_QInf << m_Beta << m_Phi << m_Ry <<m_Ctrl;

        ar<<int(m_SailAngle.size());
        for(uint is=0;is<m_SailAngle.size(); is++)
        {
            ar<< m_SailAngle[is];
        }

        m_AF.serializeFl5(ar, bIsStoring);

        ar << int(m_SailForceFF.size());
        for(uint i=0; i<m_SailForceFF.size(); i++)
        {
            ar << m_SailForceFF.at(i).xf()  << m_SailForceFF.at(i).yf()  << m_SailForceFF.at(i).zf();
            ar << m_SailForceSum.at(i).xf() << m_SailForceSum.at(i).yf() << m_SailForceSum.at(i).zf();
        }

        int N=0;
        if(isQuadMethod())           N = m_nPanel4;
        else if(isTriangleMethod())  N = 3*m_nPanel3;
        ar << N;
        for (int p=0; p<N; p++) ar << float(m_Cp[p]) << float(m_gamma[p]) << float(m_sigma[p]);


        ar << int(m_Vorton.size());
        for(uint ir=0; ir<m_Vorton.size(); ir++)
        {
            ar <<int(m_Vorton.at(ir).size());
            for(uint ic=0; ic<m_Vorton.at(ir).size(); ic++)
            {
                m_Vorton[ir][ic].serializeFl5(ar, bIsStoring);
            }
        }

        ar << int(m_VortexNeg.size());
        for(uint iv=0; iv<m_VortexNeg.size(); iv++)
        {
            m_VortexNeg[iv].serializeFl5(ar, bIsStoring);
        }


        // dynamic space allocation for the future storage of more data, without need to change the format
        nIntSpares=0;
        ar << nIntSpares;
        n=0;
        for (int i=0; i<nIntSpares; i++) ar << n;
        nDbleSpares=0;
        ar << nDbleSpares;
        dble=0.0;
        for (int i=0; i<nDbleSpares; i++) ar << dble;
    }
    else
    {
        ar >> ArchiveFormat;
        if(ArchiveFormat<100000 || ArchiveFormat>120000) return false;
        //read variables
        ar >> strange;   m_BoatName = strange.toStdString();;
        ar >> strange;   m_BtPolarName = strange.toStdString();;

        m_theStyle.serializeFl5(ar, bIsStoring);

        // ANALYSIS METHOD
        ar >> n;
        if     (n==1) m_AnalysisMethod=xfl::LLT;
        else if(n==2) m_AnalysisMethod=xfl::VLM1;
        else if(n==3) m_AnalysisMethod=xfl::VLM2;
        else if(n==4) m_AnalysisMethod=xfl::QUADS;
        else if(n==5) m_AnalysisMethod=xfl::TRILINEAR;
        else if(n==6) m_AnalysisMethod=xfl::TRIUNIFORM;

        ar >> m_bThinSurfaces;
        ar >> m_bTrefftz;

        ar >> m_bIgnoreBodyPanels;

        ar >> m_GroundHeight >> m_QInf >> m_Beta >> m_Phi;
        if(ArchiveFormat>=100005) ar >> m_Ry;
        ar >> m_Ctrl;

        ar>>n;
        m_SailAngle.resize(n);
        for(int is=0;is<n; is++)
        {
            ar>> m_SailAngle[is];
        }

        if(ArchiveFormat<100004) m_AF.serializeFl5_b17(ar, bIsStoring);
        else
        {
            if(!m_AF.serializeFl5(ar, bIsStoring)) return false;
        }
        m_AF.setOpp(0.0, m_Beta, m_Phi, m_QInf);

        ar >> n;
        for(int i=0; i<n; i++)
        {
            ar >> f >> g >> h;
            m_SailForceFF.push_back({double(f), double(g), double (h)});
            ar >> f >> g >> h;
            m_SailForceSum.push_back({double(f), double(g), double (h)});
        }

        int N=0;
        ar >> N;

        m_Cp.resize(N);
        m_gamma.resize(N);
        m_sigma.resize(N);
        for (int p=0; p<N; p++)
        {
            ar >> f >> g >> h;
            m_Cp[p]    = double(f);
            m_gamma[p]     = double(g);
            m_sigma[p] = double(h);
        }

        m_nPanel4 = m_nPanel3 = 0;
        if(isQuadMethod())           m_nPanel4 = N;
        else if(isTriangleMethod())  m_nPanel3 = N/3;

        if(ArchiveFormat>=100003)
        {
            ar >> n;
            m_Vorton.resize(n);
            for(uint ir=0; ir<m_Vorton.size(); ir++)
            {
                ar >> n;
                m_Vorton[ir].resize(n);
                for(uint ic=0; ic<m_Vorton.at(ir).size(); ic++)
                {
                    m_Vorton[ir][ic].serializeFl5(ar, bIsStoring);
                }
            }

            ar >> n;
            m_VortexNeg.resize(n);
            for(int iv=0; iv<n; iv++)
            {
                m_VortexNeg[iv].serializeFl5(ar, bIsStoring);
            }
        }

        // space allocation
        ar >> nIntSpares;
        for (int i=0; i<nIntSpares; i++) ar >> n;
        ar >> nDbleSpares;
        for (int i=0; i<nDbleSpares; i++) ar >> dble;
    }
    return true;
}


void BoatOpp::resizeResultsArrays(int N)
{
    m_Cp.resize(N);
    m_gamma.resize(N);
    m_sigma.resize(N);
    std::fill(m_Cp.begin(), m_Cp.end(), 0);
    std::fill(m_gamma.begin(), m_gamma.end(), 0);
    std::fill(m_sigma.begin(), m_sigma.end(), 0);
}


void BoatOpp::getProperties(Boat const *pBoat, double density, std::string &BOppProperties, bool bLongOutput) const
{
    std::string strong, lenunit, areaunit, forceunit, momentunit, speedunit;
    lenunit    = " " + Units::lengthUnitLabel();
    speedunit  = " " + Units::speedUnitLabel();
    forceunit  = " " + Units::forceUnitLabel();
    momentunit = " " + Units::momentUnitLabel();
    areaunit   = " " + Units::areaUnitLabel();

    double q = 0.5*density*m_QInf*m_QInf;

    BOppProperties.clear();

    if(bLongOutput)
    {
        BOppProperties += "Reference dimensions:\n";
        BOppProperties += "  area  = " + std::format("{0:9.5g}", m_AF.refArea()*Units::m2toUnit()) + areaunit + "\n";
        BOppProperties += "  chord = " + std::format("{0:9.5g}", m_AF.refChord()*Units::mtoUnit()) + lenunit+ "\n";
    }

    strong = "Ctrl    = " + std::format("{0:7.3f}", m_Ctrl);
    BOppProperties += strong + "\n";

    strong = "AWS" + INFch + "    = " + std::format("{0:7.3f}", m_QInf);
    BOppProperties += strong + speedunit+"\n";

    strong = "AWA" + INFch + "    = " + std::format("{0:7.3f}", m_Beta);
    BOppProperties += strong +DEGch+"\n";

    strong = PHIch + "       = " + std::format("{0:7.3f}", m_Phi);
    BOppProperties += strong +DEGch+"\n";

    for(uint is=0; is<m_SailAngle.size(); is++)
    {
        strong = std::format("Sail_{0:d} Angle = ", is+1) + std::format("{0:9.3g}", m_SailAngle.at(is));
        BOppProperties += strong + DEGch+"\n";
    }

    BOppProperties += "CL           = " + std::format("{0:9.5g}", m_AF.CSide())+"\n";
    BOppProperties += "CD           = " + std::format("{0:9.5g}", m_AF.CD())+"\n";
    if(bLongOutput)
    {
        BOppProperties += "   CDi       = " + std::format("{0:9.5g}", m_AF.CDi())+"\n";
        BOppProperties += "   CDv       = " + std::format("{0:9.5g}", m_AF.CDv())+"\n";
    }
    BOppProperties += "Cx           = " + std::format("{0:9.5g}", m_AF.Cx())+"\n";
    BOppProperties += "Cy           = " + std::format("{0:9.5g}", m_AF.Cy())+"\n";

    strong = "Far Field Fx = " + std::format("{0:9.5g}", m_AF.fffx() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";
    strong = "Far Field Fy = " + std::format("{0:9.5g}", m_AF.fffy() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";
    strong = "Far Field Fz = " + std::format("{0:9.5g}", m_AF.fffz() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";

    strong = "Summed Fx    = " + std::format("{0:9.5g}", m_AF.fsumx() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";
    strong = "Summed Fy    = " + std::format("{0:9.5g}", m_AF.fsumy() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";
    strong = "Summed Fz    = " + std::format("{0:9.5g}", m_AF.fsumz() * q * Units::NtoUnit());
    BOppProperties += strong +forceunit+"\n";

    strong  = "Mx = " + std::format("{0:9.5g}", (m_AF.Mi()+m_AF.Mv()).x * q * Units::NmtoUnit());
    BOppProperties += strong +momentunit+"\n";
    strong  = "My = " + std::format("{0:9.5g}", (m_AF.Mi()+m_AF.Mv()).y * q * Units::NmtoUnit());
    BOppProperties += strong +momentunit+"\n";
    strong  = "Mz = " + std::format("{0:9.5g}", (m_AF.Mi()+m_AF.Mv()).z * q * Units::NmtoUnit());
    BOppProperties += strong +momentunit;

    if(!bLongOutput) return;
    BOppProperties +="\n";

    strong  = "XCE = Mz/Fy =" + std::format("{0:9.5g}", (m_AF.Mi()+m_AF.Mv()).z /m_AF.fsumy() * Units::mtoUnit());
    BOppProperties += strong +lenunit+"\n";
    strong  = "ZCE = Mx/Fy = " + std::format("{0:9.5g}", (m_AF.Mi()+m_AF.Mv()).x /m_AF.fsumy() * Units::mtoUnit());
    BOppProperties += strong +lenunit +"\n";

    for(uint is=0; is<m_SailForceFF.size(); is++)
    {
        BOppProperties += pBoat->sailAt(is)->name() + ": \n";
        strong = "   Far Field Fx= " + std::format("{0:9.5g}", m_SailForceFF.at(is).x * q * Units::NtoUnit());
        BOppProperties += strong +forceunit+"\n";
        strong = "   Far Field Fy= " + std::format("{0:9.5g}", m_SailForceFF.at(is).y * q * Units::NtoUnit());
        BOppProperties += strong +forceunit+"\n";
        strong = "   Far Field Fz= " + std::format("{0:9.5g}", m_SailForceFF.at(is).z * q * Units::NtoUnit());
        BOppProperties += strong +forceunit+"\n";
        strong = "   Summed Fx   = " + std::format("{0:9.5g}", m_SailForceSum.at(is).x * q * Units::NtoUnit());
        BOppProperties += strong +forceunit+"\n";
        strong = "   Summed Fy   = " + std::format("{0:9.5g}", m_SailForceSum.at(is).y * q * Units::NtoUnit());
        BOppProperties += strong +forceunit+"\n";
        strong = "   Summed Fz   = " + std::format("{0:9.5g}", m_SailForceSum.at(is).z * q * Units::NtoUnit());
        BOppProperties += strong +forceunit +"\n";
    }
    BOppProperties.pop_back(); //last \n
}


Vector3d BoatOpp::windDir() const
{
    return objects::windDirection(0.0, -m_Beta);
}


std::string BoatOpp::title(bool bLong) const
{
    std::string strong;

    if(bLong)
    {
        if(!isLLTMethod()) //always
        {
            if(isVLM1())
            {
                strong += "VLM1";
            }
            else if(isVLM2())
            {
                strong += "VLM2";
            }
            else if(isQuadMethod())
            {
                strong += "Quads";
            }
            else if(isTriUniformMethod())
            {
                strong += "TriUni";
            }
            else if(isTriLinearMethod())
            {
                strong += "TriLin";
            }
        }

        strong +="-";
    }

    strong += std::format("{0:5.2f}-", ctrl());

    if(fabs(beta())>PRECISION) strong += std::format("{0:5.2f}°-", beta());

    return strong;
}



void BoatOpp::exportMainDataToString(Boat const*, std::string &poppdata, xfl::enumTextFileType filetype, std::string const &textsep) const
{
    std::string strange;
    std::string title;
    std::string len = Units::lengthUnitLabel();
    std::string inertia = Units::inertiaUnitLabel();

    std::string sep = "  ";
    if(filetype==xfl::CSV) sep = textsep+ " ";

    poppdata += boatName()+"\n";
    poppdata += polarName()+"\n\n";
    poppdata += "ctrl       " + sep
              + "beta       " + sep
              + "phi        " + sep
              + "VInf("+Units::speedUnitLabel()+")" +sep
              + "h("+Units::lengthUnitLabel()+")\n";
    strange = std::format("{:11.5g}", m_Ctrl);
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_Beta);
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_Phi);
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_QInf*Units::mstoUnit());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_GroundHeight*Units::mtoUnit());
    poppdata += strange;

    poppdata += "\n\n";

    strange = "CL         " + sep
            + "CX         " + sep
            + "CY         " + sep
            + "CD_inviscid" + sep
            + "CD_viscous " + sep
            + "Cl         " + sep
            + "Cm_inviscid" + sep
            + "Cm_viscous " + sep
            + "Cn_inviscid" + sep
            + "Cn_viscous\n";
    poppdata += strange;

    strange = std::format("{:11.5g}", m_AF.CL());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cx());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cy());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.CDi());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.CDv());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cli());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cmi());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cmv());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cni());
    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.Cnv());
    poppdata += strange;

    poppdata += "\n";

    strange = "CP.x("+len+")";      for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CP.y("+len+")";      for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CP.z("+len+")";      for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "NP.x("+len+")";      for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    poppdata += "\n";


    strange = std::format("{:11.5g}", m_AF.centreOfPressure().x*Units::mtoUnit());    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.centreOfPressure().y*Units::mtoUnit());    poppdata += strange+sep;
    strange = std::format("{:11.5g}", m_AF.centreOfPressure().z*Units::mtoUnit());    poppdata += strange+sep;
    poppdata += strange + "\n\n";

    strange = "mass("+Units::massUnitLabel()+")";  for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG.x("+len+")";                    for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG.y("+len+")";                    for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG.z("+len+")";                    for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG_Ixx("+inertia+")";              for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG_Iyy("+inertia+")";              for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG_Izz("+inertia+")";              for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    strange = "CoG_Ixz("+inertia+")";              for(int i=int(strange.length()); i<11; i++) strange+=" ";   poppdata += strange+sep;
    poppdata +="\n";


    poppdata += strange + "\n\n";


    title = "y("+Units::lengthUnitLabel()+")";  for(int i=int(title.length()); i<11; i++) title+=" ";   title +=sep;
    title += "Re         " + sep;
    title += "Ai         " + sep;
    title += "Cd_i       " + sep;
    title += "Cd_v       " + sep;
    title += "Cl         " + sep;
    title += "CP.x(%)    " + sep;
    title += "Trans.top  " + sep;
    title += "Trans.bot  " + sep;
    title += "Cm_i       " + sep;
    title += "Cm_v       " + sep;
    title += "Bending.mom" + sep;
    title += "Vd.x       " + sep;
    title += "Vd.y       " + sep;
    title += "Vd.z       " + sep;
    title += "F.x        " + sep;
    title += "F.y        " + sep;
    title += "F.z        ";
    title += '\n';
}


void BoatOpp::exportPanel3DataToString(Boat const*pBoat,
                                       xfl::enumTextFileType exporttype,
                                       std::string &paneldata) const
{
    std::string strong;

    if(exporttype==xfl::TXT) paneldata += " Panel        CtrlPt.x        CtrlPt.y        CtrlPt.z           Nx               Ny             Nz            Area             Cp\n";
    else                     paneldata += "Panel,CtrlPt.x,CtrlPt.y,CtrlPt.z,Nx,Ny,Nz,Area,Cp\n";


    for(int iw=0; iw<pBoat->nSails(); iw++)
    {
        Sail const *pSail = pBoat->sailAt(iw);

        paneldata += pSail->name()+ "_Cp Coefficients"+"\n";
        int p=pSail->m_FirstPanel3Index;

        for(int i3=0; i3<pSail->nPanel3(); i3++)
        {
            Panel3 const &p3 = pBoat->triMesh().panelAt(p);

            double cp=0;
            for(int in=0; in<3; in++) cp += m_Cp.at(p3.index()*3+in);
            cp /= 3.0;

            if(exporttype==xfl::TXT)
                strong = std::format("{:d}     {:11g}     {:11g}     {:11g}     {:11g}     {:11g}     {:11g}     {:11g}     {:11g}\n",
                                     p, p3.CoG().x,  p3.CoG().y, p3.CoG().z,  p3.normal().x, p3.normal().y, p3.normal().z, p3.area(), cp);
            else
                strong = std::format("{:d}, {:11g}, {:11g}, {:11g}, {:11g}, {:11g}, {:11g}, {:11g}, {:11g}\n",
                                     p, p3.CoG().x,  p3.CoG().y, p3.CoG().z,  p3.normal().x, p3.normal().y, p3.normal().z, p3.area(), cp);

            paneldata += strong;
            p++;
        }

        paneldata += ("\n\n");
    }
}
