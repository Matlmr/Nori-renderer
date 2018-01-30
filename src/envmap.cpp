/*
 This file is part of Nori, a simple educational ray tracer
 
 Copyright (c) 2015 by Romain Prévost
 
 Nori is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License Version 3
 as published by the Free Software Foundation.
 
 Nori is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <nori/emitter.h>
#include <nori/bitmap.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

// Same typedef than in Bitmap
typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix;

class EnvironmentMap : public Emitter {
public:
    EnvironmentMap(const PropertyList &props) {
        std::string filename = props.getString("mapname");
        m_bitmap = Bitmap(filename);
        m_cols = m_bitmap.cols();
        m_rows = m_bitmap.rows();
        preprocessing();
    }
    
    virtual bool isEnvEmitter() const {
        return true;
    }
    
    virtual std::string toString() const override {
        return tfm::format(
                           "EnvironmentMap[\n"
                           );
    }
    
    virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
        
        Point2f uv = sphericalCoordinates(lRec.wi.normalized());
        
        float row = uv.x() * INV_PI * (m_rows - 1);
        float col = uv.y() * 0.5f * INV_PI * (m_cols - 1);
        
        return bilinearInterpolation(row, col);
    }
    
    virtual float pdf(const EmitterQueryRecord &lRec) const override {
        Point2f uv = sphericalCoordinates(lRec.wi.normalized());
        float row = uv.x() * INV_PI * (m_rows - 1);
        float col = uv.y() * 0.5f * INV_PI * (m_cols - 1);
        
        int x = (int) row;
        int y = (int) col;
        if (x > m_rows - 1) x = m_rows - 1;
        if (y > m_cols - 1) y = m_cols - 1;
        if (x < 0) x = 0;
        if (y < 0) y = 0;

        return m_mPDF(0,x) * m_cPDF(x, y);
        return 1.f;
    }
    
    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
        
        float u, v, pdf_u, pdf_v;
        sample1D(0, m_mPDF, m_mCDF, sample.x(), u, pdf_u);
        sample1D(u, m_cPDF, m_cCDF, sample.y(), v, pdf_v);
        u *= M_PI / (m_rows - 1);
        v *= 2 * M_PI / (m_cols - 1);
        Vector3f w = Vector3f(sin(u) * cos(v), sin(u) * sin(v), cos(u)).normalized();
        
        lRec.wi = w;
        float detJac = (m_rows - 1) * (m_cols - 1) * 0.5f * pow(INV_PI,2) / Frame::sinTheta(lRec.wi);
        
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, 100000);
        
        float pdfv = pdf(lRec) * detJac;
        //return eval(lRec);
        //return eval(lRec) / pdfv;
        return eval(lRec) / pdf(lRec) / detJac;
        
        
    }
    
    void preprocessing () {
        m_cPDF = matrix(m_rows,m_cols);
        m_cCDF = matrix(m_rows,m_cols + 1);
        m_mPDF = matrix(1,m_rows);
        m_mCDF = matrix(1,m_rows + 1);
        m_luminance = matrix(m_rows, m_cols);
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_luminance(i,j) = sqrt(0.299*pow(m_bitmap(i,j).r(),1) + 0.587*pow(m_bitmap(i,j).g(),1) + 0.114*pow(m_bitmap(i,j).b(),1)) + Epsilon / 1000000;
            }
        }
        matrix sum(1, m_rows);
        for (int i = 0; i < m_rows; ++i) precompute1D(i, m_luminance, m_cPDF, m_cCDF, sum(0,i));
        float I;
        precompute1D(0, sum, m_mPDF, m_mCDF, I);
    }
    
    void precompute1D(int row, matrix &f, matrix &pf, matrix &Pf, float &I) {
        int nf = f.cols();
        I = 0.f;
        for (int i = 0; i < nf; ++i) I += f(row,i);
        if (I == 0.0f) return;
        for (int i = 0; i < nf; ++i) pf(row,i) = f(row,i) / I;
        Pf(row,0) = 0.f;
        for (int i = 1; i < nf; ++i) Pf(row,i) = Pf(row,i-1) + pf(row,i - 1);
        Pf(row,nf) = 1.f;
    }
    
    void sample1D(int row, const matrix &pf, const matrix &Pf, float unif, float &x, float &pdf) const {
        int i;
        for (i = 0; i < Pf.cols(); ++i) {
            if (unif >= Pf(row, i) && unif < Pf(row, i+1)) break;
        }
        
        float t = (Pf(row,i + 1) - unif) / (Pf(row,i+1) - Pf(row,i));
        x = (1 - t) * i + t * (i + 1);
        pdf = pf(row,i);
    }
    
    Color3f bilinearInterpolation(float x, float y) const {
        
        int x1 = (int) x;
        int y1 = (int) y;
        int x2 = x1 + 1;
        int y2 = y1 + 1;
        Color3f Q11 = 0.f;
        if (x1 >= 0 && x1 < m_rows && y1 >= 0 && y1 < m_cols) Q11 = m_bitmap(x1, y1);
        Color3f Q12 = 0.f;
        if (x1 >= 0 && x1 < m_rows && y2 >= 0 && y2 < m_cols) Q12 = m_bitmap(x1, y2);
        Color3f Q21 = 0.f;
        if (x2 >= 0 && x2 < m_rows && y1 >= 0 && y1 < m_cols) Q21 = m_bitmap(x2, y1);
        Color3f Q22 = 0.f;
        if (x2 >= 0 && x2 < m_rows && y2 >= 0 && y2 < m_cols) Q22 = m_bitmap(x2, y2);
        int Dx = x2 - x1;
        int Dy = y2 - y1;
        float dx2 = x2 - x;
        float dy2 = y2 - y;
        float dy1 = y - y1;
        float dx1 = x - x1;
        if (Dx == 0.f || Dy == 0.f) return 0.f;
        return 1.0 / (Dx * Dy) * (Q11 * dx2 * dy2 +
                                  Q21 * dx1 * dy2 +
                                  Q12 * dx2 * dy1 +
                                  Q22 * dx1 * dy1);
    }
    
protected:
    Bitmap m_bitmap;
    int m_cols, m_rows;
    matrix m_luminance;
    matrix m_mPDF, m_mCDF, m_cPDF, m_cCDF;
};

NORI_REGISTER_CLASS(EnvironmentMap, "envmap")
NORI_NAMESPACE_END
