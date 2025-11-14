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


#pragma once


#include <condition_variable>
#include <mutex>
#include <queue>



#include <api/fl5lib_global.h>


#include <xfoil.h>
#include <api/analysisrange.h>
#include <api/utils.h>

class Foil;
class Polar;
class OpPoint;


struct FoilAnalysis
{
    Foil *pFoil{nullptr};            /**< a pointer to the Foil object to be analyzed by the thread */
    Polar *pPolar{nullptr};          /**< a pointer to the Polar object to be analyzed by the thread */
    std::vector<AnalysisRange> range;
};


/**
* @class  XFoilTask
* This file implements the management task of an XFoil calculation.
*/
class FL5LIB_EXPORT XFoilTask
{
    public:
        XFoilTask();


        void run();

        void setAnalysisStatus(xfl::enumAnalysisStatus status) {m_AnalysisStatus = status;}
        bool isCancelled() const {return m_AnalysisStatus==xfl::CANCELLED;}
        bool isRunning()   const {return m_AnalysisStatus==xfl::RUNNING;}
        bool isPending()   const {return m_AnalysisStatus==xfl::PENDING;}
        bool isFinished()  const {return m_AnalysisStatus==xfl::FINISHED || m_AnalysisStatus==xfl::CANCELLED;}

        bool hasErrors()   const {return m_bErrors;}

        Foil const *foil() const {return m_pFoil;}
        Polar const *polar() const {return m_pPolar;}

        std::vector<OpPoint*> const &operatingPoints() const {return m_OpPoints;}

        bool processCl(double Cl, double Re, double &Cd, double &XTrTop, double &XTrBot, bool &bCv);
        bool processClList(const std::vector<double> &ClList, const std::vector<double> &ReList, std::vector<double> &CdList, std::vector<double> &XTrTopList, std::vector<double> &XTrBotList, std::vector<bool> &CvList);

        void clearRanges() {m_AnalysisRange.clear();}
        void setAnalysisRanges(std::vector<AnalysisRange> const &ranges) {m_AnalysisRange=ranges;}
        void appendRange(AnalysisRange const &range) {m_AnalysisRange.push_back(range);}
        void appendRange(bool b, double vmin, double vmax, double vinc) {m_AnalysisRange.push_back({b, vmin, vmax, vinc});}

        bool initialize(FoilAnalysis *pFoilAnalysis, bool bStoreOpp, bool bViscous=true, bool bInitBL=true);
        bool initialize(Foil *pFoil, Polar *pPolar, bool bStoreOpp, bool bViscous=true, bool bInitBL=true);

        XFoil const &XFoilInstance() const {return m_XFoilInstance;}

        void setAlphaRange(double vMin, double vMax, double vDelta);
        void setClRange(double vMin, double vMax, double vDelta);


        void traceLog(const std::string &str);

        static void cancelAnalyses() {s_bCancel=true;}
        static void setCancelled(bool b) {s_bCancel=b;}

        static void setSkipOpp(bool b) {s_bSkipOpp=b;}
        static bool bSkipOpp() {return s_bSkipOpp;}

        static void setSkipPolar(bool b) {s_bSkipPolar=b;}
        static bool bSkipPolar() {return s_bSkipPolar;}

        static int maxIterations() {return s_IterLim;}
        static void setMaxIterations(int maxiter) {s_IterLim=maxiter;}

        static bool bAutoInitBL() {return s_bAutoInitBL;}
        static void setAutoInitBL(bool b) {s_bAutoInitBL=b;}

        static bool bStoreOpps() {return s_bStoreOpp;}
        static void setStoreOpps(bool b) {s_bStoreOpp=b;}

        static bool bAlpha()     {return s_bAlpha;}
        static void setAoAAnalysis(bool b) {s_bAlpha=b;}

        static bool bViscous()   {return s_bViscous;}
        static void setViscous(bool b) {s_bViscous=b;}

        static bool bInitBL()    {return s_bInitBL;}
        static void setInitBL(bool b) {s_bInitBL=b;}

        static double CdError() {return s_CdError;}
        static double setCdError(double cderr) {return s_CdError=cderr;}

    protected:
        int loop();
        bool alphaSequence(bool bAlpha);
        bool thetaSequence();
        bool ReSequence();
        void addXFoilData(OpPoint *pOpp, XFoil &xfoil, const Foil *pFoil);


    private:
        bool m_bErrors;


        XFoil m_XFoilInstance;     /**< An instance of the XFoil class specific for this object */

        std::vector<OpPoint*> m_OpPoints;

        Foil *m_pFoil;                 /**< A pointer to the instance of the Foil object for which the calculation is performed */
        Polar *m_pPolar;                /**< A pointer to the instance of the Polar object for which the calculation is performed */

        xfl::enumAnalysisStatus m_AnalysisStatus;

        std::vector<AnalysisRange> m_AnalysisRange;

        static bool s_bSkipPolar;
        static bool s_bCancel;            /**< True if the user has asked to cancel the analysis */
        static bool s_bSkipOpp;

        static int  s_IterLim;
        static bool s_bAutoInitBL;        /**< true if the BL initialization is left to the code's decision */
        static bool s_bViscous;           /**< true if performing a viscous calculation, false if inviscid */
        static bool s_bAlpha;             /**< true if performing an analysis based on aoa, false if based on Cl */
        static bool s_bInitBL;            /**< true if the boundary layer should be initialized for the next xfoil calculation */
        static bool s_bStoreOpp;          /**< true if the operating points should be stored */
        static double s_CdError;          /**< discard points with |Cd| less than this value: these operating points are likely erroneous (spurious?) */


    public:
        // thread related variables to share the message queue with the calling thread
        std::mutex m_mtx;
        std::condition_variable m_cv;
        std::queue<std::string> m_theMsgQueue;
};




