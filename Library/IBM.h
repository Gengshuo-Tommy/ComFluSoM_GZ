/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2021 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License fo`r more details.                        *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// IBM methods

#ifndef IBM_H
#define IBM_H

struct StructMeshIndexWeight
{
    Vector3i Index;
    double W = 0.;
};

void CalWeightStructMesh(size_t dim, int domSize[3], bool periodic[3], size_t kernel, Vector3d& x, vector<StructMeshIndexWeight>& l)
{
    l.resize(0);

    double boxL = 0.;
    double (*DeltaFunc)(double s) = NULL;

    switch(kernel)
    {
        case 0:
        {
            boxL = 1.51;
            DeltaFunc = &Delta3P;
        }break;
        case 1:
        {
            boxL = 2.01;
            DeltaFunc = &Delta4P;
        }break;
    }
    // Find node
    Vector3i min0 (0,0,0);
    Vector3i max0 (0,0,0);

    for(size_t d=0; d<dim; ++d)
    {
        min0(d) = floor(x(d)-boxL);
        max0(d) = floor(x(d)+boxL);
    }

    double sumw = 0.;

    for (int k=min0(2); k<=max0(2); ++k)
    for (int j=min0(1); j<=max0(1); ++j)
    for (int i=min0(0); i<=max0(0); ++i)
    {
        Vector3d xn (i,j,k);
        double s = (x-xn).norm();
        double w = DeltaFunc(s);
        if (w>0.)
        {
            bool withInDomain = true;
            Vector3i xni (i,j,k);
            for(size_t d=0; d<dim; ++d)
            {
                if (periodic[d])    xni[d] = (xni[d]+domSize[d]+1)%(domSize[d]+1);
                if (xni[d]<0 || xni[d]>domSize[d])
                {
                    withInDomain = false;
                    break;
                }
            }

            if (withInDomain)
            {
                sumw += w;
                StructMeshIndexWeight indw;
                indw.Index = xni;
                indw.W = w;
                l.push_back(indw);
            }
        }
    }

    for (size_t m=0; m<l.size(); ++m)   l[m].W /= sumw;
}

// T1 need support += /=, the user need to intialize "value" with zero by themself
template <typename T0, typename T1>
void InterpolateStructMesh(vector<StructMeshIndexWeight>& l, T1& meshValue, T0& value)
{
    for (size_t m=0; m<l.size(); ++m)
    {
        size_t i = l[m].Index(0);
        size_t j = l[m].Index(1);
        size_t k = l[m].Index(2);
        double w = l[m].W;     
        value += w*meshValue[i][j][k];   
    }
}

// need to parallelize this part with lock
template <typename T0, typename T1>
void DistributeStructMesh(vector<StructMeshIndexWeight>& l, T0& value, T1& meshValue)
{
    for (size_t m=0; m<l.size(); ++m)
    {
        size_t i = l[m].Index(0);
        size_t j = l[m].Index(1);
        size_t k = l[m].Index(2);
        double w = l[m].W;

        meshValue[i][j][k] += w*value;
    }
}

// distribute two variables at same time for efficiency
// need to parallelize this part with lock
template <typename T0, typename T1, typename T2, typename T3>
void DistributeStructMesh(vector<StructMeshIndexWeight>& l, T0& value0, T1& meshValue0, T2& value1, T3& meshValue1)
{
    for (size_t m=0; m<l.size(); ++m)
    {
        size_t i = l[m].Index(0);
        size_t j = l[m].Index(1);
        size_t k = l[m].Index(2);
        double w = l[m].W;

        meshValue0[i][j][k] += w*value0;
        meshValue1[i][j][k] += w*value1;
    }
}

#endif
