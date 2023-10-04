#pragma once

#include "animaSPANClusteringFilter.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
//#include <itkSymmetricEigenAnalysis.h>

namespace anima {
	
	template <class DataType, unsigned int DataDimension>
    SPANClusteringFilter <DataType, DataDimension>::
    SPANClusteringFilter()
	{
		m_ClassMemberships.clear();
		m_ReverseClassMemberships.clear();
		m_InputData.clear();
		m_ConnectivityMatrix.clear();
		m_ModularityMatrix.clear();
		m_DegreeMatrix.clear();
		m_PreventClusterFromSplitting.clear();
		
		m_OrientationData = false;

		m_NbInputs = 0;
		m_NumberOfClusters = 0;
		
		m_TotalNumberOfEdges = 0.0;
	}
	
	template <class DataType, unsigned int DataDimension>
	SPANClusteringFilter <DataType, DataDimension>::
    ~SPANClusteringFilter()
	{
	}
	
    template <class DataType, unsigned int DataDimension>
    void
    SPANClusteringFilter <DataType, DataDimension>::
    SetInputData(DataHolderType &data)
	{
		if (data.size() == 0)
			return;
		
		m_InputData = data;
		
		m_NbInputs = m_InputData.size();
	}
	
	template <class DataType, unsigned int DataDimension>
	void
    SPANClusteringFilter <DataType,DataDimension>::
    InitializeClassMemberships()
	{
		m_NumberOfClusters = 1;
		
		m_ClassMemberships.resize(m_NbInputs);
		m_ReverseClassMemberships.resize(m_NumberOfClusters);
		m_ReverseClassMemberships[0].resize(m_NbInputs);
        for (unsigned int i = 0;i < m_NbInputs;++i)
        {
        	m_ClassMemberships[i] = 0;
        	m_ReverseClassMemberships[0][i] = i;
        }
        
        m_PreventClusterFromSplitting.resize(m_NumberOfClusters);
        m_PreventClusterFromSplitting[0] = false;
	}
	
	template <class DataType, unsigned int DataDimension>
	typename SPANClusteringFilter <DataType, DataDimension>::MatrixType
	SPANClusteringFilter <DataType, DataDimension>::
	ComputeCMatrix(unsigned int parentClusterIndex, VectorType &BtPerOne)
	{
		unsigned int parentClusterSize = m_ReverseClassMemberships[parentClusterIndex].size();

		MatrixType C(parentClusterSize,DataDimension);
		C.fill(0.0);
		
		BtPerOne.set_size(DataDimension);
		BtPerOne.fill(0);
		
		for (unsigned int j = 0;j < DataDimension;++j)
		{
			BtPerOne(j) = 0;
			for (unsigned int i = 0;i < parentClusterSize;++i)
			{
				unsigned int currentParticle = m_ReverseClassMemberships[parentClusterIndex][i];
				BtPerOne(j) += m_InputData[currentParticle][j];
			}
		}
		
		for (unsigned int i = 0;i < parentClusterSize;++i)
		{
			unsigned int currentParticle = m_ReverseClassMemberships[parentClusterIndex][i];
			
			double tmpVal = 0;
			
			for (unsigned int j = 0;j < DataDimension;++j)
				tmpVal += m_InputData[currentParticle][j] * BtPerOne(j);
				
			if (tmpVal <= 0)
			{
				std::cerr << "SPAN Clustering: negative entries in similarity matrix; this might have also happened before this error displays." << std::endl;
				exit(-1);
			}
			
			for (unsigned int j = 0;j < DataDimension;++j)
				C(i,j) = m_InputData[currentParticle][j] / std::sqrt(tmpVal);
		}
		
		return C;
	}
	
	template <class DataType, unsigned int DataDimension>
	void
	SPANClusteringFilter <DataType, DataDimension>::
	FastBipartioning(MatrixType &data, unsigned int parentIndex, MembershipType &members1, MembershipType &members2)
	{	
		unsigned int parentClusterSize = data.rows();
		vnl_svd < double > svd(data);
		
		unsigned int pos = 0;
		while (svd.W(pos) == svd.W(0))
			++pos;
		
		VectorType tmpVec = svd.U().get_column(pos);
		
		for (unsigned int i = 0;i < parentClusterSize;++i)
		{
			unsigned int currentParticle = m_ReverseClassMemberships[parentIndex][i];
					
			if (tmpVec(i) < 0)
				members1.push_back(currentParticle);
			else
				members2.push_back(currentParticle);
		}
	}
	
	template <class DataType, unsigned int DataDimension>
	double
	SPANClusteringFilter <DataType, DataDimension>::
	GetNewmanGirvanModularity(MembershipType &members1, MembershipType &members2, VectorType &BtPerOne)
	{
		unsigned int numGrp1 = members1.size();
		unsigned int numGrp2 = members2.size();
		unsigned int numData = numGrp1 + numGrp2;
		
		double L = -(double)numData, L1 = -(double)numGrp1, L2 = -(double)numGrp2;
		double O11 = -(double)numGrp1, O22 = -(double)numGrp2;
		
		for (unsigned int i = 0;i < DataDimension;++i)
		{
			L += BtPerOne(i) * BtPerOne(i);
		
			double tmpVal1 = 0;
			for (unsigned int j = 0;j < numGrp1;++j)
			{
				unsigned int currentParticle = members1[j];
				tmpVal1 += m_InputData[currentParticle][i];
			}
		
			double tmpVal2 = 0;
			for (unsigned int j = 0;j < numGrp2;++j)
			{
				unsigned int currentParticle = members2[j];
				tmpVal2 += m_InputData[currentParticle][i];
			}
		
			L1 += tmpVal1 * BtPerOne(i);
			L2 += tmpVal2 * BtPerOne(i);
			O11 += tmpVal1 * tmpVal1;
			O22 += tmpVal2 * tmpVal2;
		}
		
		L1 /= L;
		L2 /= L;
		O11 /= L;
		O22 /= L;
		
		double resVal = O11 + O22 - (L1 * L1 + L2 * L2);
		return resVal; 
	}
	
	template <class DataType, unsigned int DataDimension>
	void
	SPANClusteringFilter <DataType, DataDimension>::
	SplitCluster(unsigned int parentCluster)
	{
		unsigned int parentClusterSize = m_ReverseClassMemberships[parentCluster].size();
		
//		if (parentClusterSize == 1)
//		{
//			m_PreventClusterFromSplitting[parentCluster] = true;
//			return;
//		}
		
		MembershipType members1, members2;
		double Q = 0;
		
		if (m_OrientationData)
		{
			// perform traditional spectral clustering
			Q = this->SpectralClustering(parentCluster, members1, members2);
		}
		else
		{
			// Get C matrix for the subset B matrix stored in data
			VectorType BtPerOne(DataDimension,0);
			MatrixType C = this->ComputeCMatrix(parentCluster, BtPerOne);
			
			// Split data in two according to fast bipartioning method
			this->FastBipartioning(C, parentCluster, members1, members2);
			
			// Compute Newman-Girvan modularity to test validity of the split
			Q = this->GetNewmanGirvanModularity(members1, members2, BtPerOne);
		}
		
//		if (Q < -0.5 || Q > 1)
//		{
//			std::cerr << "Modularity out of bounds: " << Q << std::endl;
//			exit(-1);
//		}

		if (Q > 1.0e-5 && members1.size() != parentClusterSize && members2.size() != parentClusterSize)
		{
			// Update class memberships		
			for (unsigned int i = 0;i < members2.size();++i)
				m_ClassMemberships[members2[i]] = m_NumberOfClusters;
			
			m_ReverseClassMemberships[parentCluster] = members1;
			m_ReverseClassMemberships.push_back(members2);
			++m_NumberOfClusters;
			
			m_PreventClusterFromSplitting.push_back(false);
		}
		else
			m_PreventClusterFromSplitting[parentCluster] = true;
	}
	
	template <class DataType, unsigned int DataDimension>
	void
    SPANClusteringFilter <DataType,DataDimension>::
    Update()
	{
		if (m_NbInputs == 0)
		{
			std::cerr << "Input data is empty." << std::endl;
			exit(-1);
		}
			
		this->InitializeClassMemberships();
		
		if (m_OrientationData)
		{
			if (m_NbInputs == 2)
			{
				double t = 0;
				for (unsigned int i = 0;i < DataDimension;++i)
					t += m_InputData[0][i] * m_InputData[1][i];
					
				if (t * t < 0.5)
				{
					++m_NumberOfClusters;
					m_ClassMemberships[1] = 1;
					m_ReverseClassMemberships[0].pop_back();
					m_ReverseClassMemberships.push_back(m_ReverseClassMemberships[0]);
					m_ReverseClassMemberships[1][0] = 1;
				}
				
				return;
			}
			this->GetConnectivityMatrix();
			this->GetModularityMatrix();
		}
		
		unsigned int oldNumberOfClusters = 0;
		while (m_NumberOfClusters > oldNumberOfClusters)
		{
			oldNumberOfClusters = m_NumberOfClusters;
			
			for (unsigned int i = 0;i < oldNumberOfClusters;++i)
				if (!m_PreventClusterFromSplitting[i])
					this->SplitCluster(i);
		}
	}
	
	template <class DataType, unsigned int DataDimension>
	void
    SPANClusteringFilter <DataType,DataDimension>::
    GetConnectivityMatrix()
	{
		m_ConnectivityMatrix.set_size(m_NbInputs,m_NbInputs);
		m_ConnectivityMatrix.fill(0.0);
		
		m_DegreeMatrix.set_size(m_NbInputs);
		m_DegreeMatrix.fill(0.0);
		
		m_TotalNumberOfEdges = 0.0;
		
		for (unsigned int i = 0;i < m_NbInputs;++i)
		{
			for (unsigned int j = i+1;j < m_NbInputs;++j)
			{
				double t = 0;
				for (unsigned int k = 0;k < DataDimension;++k)
					t += m_InputData[i][k] * m_InputData[j][k];
				
				t *= t;
					
				m_ConnectivityMatrix(i,j) = t;
				m_ConnectivityMatrix(j,i) = t;
				m_DegreeMatrix(i) += t;
				m_TotalNumberOfEdges += t;	
			}
			
			for (unsigned int j = 0;j < i;++j)
				m_DegreeMatrix(i) += m_ConnectivityMatrix(j,i);	
		}
	}
	
	template <class DataType, unsigned int DataDimension>
	void
    SPANClusteringFilter <DataType,DataDimension>::
    GetModularityMatrix()
	{
		m_ModularityMatrix.set_size(m_NbInputs,m_NbInputs);
		
		for (unsigned int i = 0;i < m_NbInputs;++i)
		{
			for (unsigned int j = i;j < m_NbInputs;++j)
			{
				double t = m_ConnectivityMatrix(i,j) - m_DegreeMatrix(i) * m_DegreeMatrix(j) / (2.0 * m_TotalNumberOfEdges);
				m_ModularityMatrix(i,j) = t;
					
				if (i != j)
					m_ModularityMatrix(j,i) = t;
			}
		}
	}
	
	template <class DataType, unsigned int DataDimension>
	typename SPANClusteringFilter <DataType,DataDimension>::MatrixType
    SPANClusteringFilter <DataType,DataDimension>::
    GetSubModularityMatrix(unsigned int parentIndex)
	{
		MembershipType memberships = m_ReverseClassMemberships[parentIndex];
		unsigned int numData = memberships.size();
		
		VectorType subDegreeMatrix(numData,0.0);
		
		for (unsigned int i = 0;i < numData;++i)
			for (unsigned int j = 0;j < numData;++j)
				subDegreeMatrix(i) += m_ModularityMatrix(memberships[i],memberships[j]);
		
		MatrixType B(numData,numData);
		
		for (unsigned int i = 0;i < numData;++i)
		{
			for (unsigned int j = i;j < numData;++j)
			{
				B(i,j) = m_ModularityMatrix(memberships[i],memberships[j]);
				
				if (i == j)
					B(i,j) -= subDegreeMatrix(i);
				else
					B(j,i) = B(i,j);
			}
		}
		
		return B;
	}
	
	template <class DataType, unsigned int DataDimension>
	double
    SPANClusteringFilter <DataType,DataDimension>::
    SpectralClustering(unsigned int parentIndex, MembershipType &members1, MembershipType &members2)
	{	
		unsigned int numData = m_ReverseClassMemberships[parentIndex].size();
		
		MatrixType B = this->GetSubModularityMatrix(parentIndex);
		
		vnl_symmetric_eigensystem <double> eigSystem(B);
		VectorType tmpVec = eigSystem.get_eigenvector(numData-1);
		
		for (unsigned int i = 0;i < numData;++i)
		{
			unsigned int currentParticle = m_ReverseClassMemberships[parentIndex][i];
					
			if (tmpVec(i) < 0)
				members1.push_back(currentParticle);
			else
				members2.push_back(currentParticle);
		}
		
		// Compute modularity for splitting this cluster into two subclusters
		double s_i, s_j, Q = 0;
		for (unsigned int i = 0;i < numData;++i)
		{
			if (tmpVec(i) < 0)
				s_i = -1;
			else
				s_i = 1;
				
			for (unsigned int j = 0;j < numData;++j)
			{
				if (tmpVec(j) < 0)
					s_j = -1;
				else
					s_j = 1;
					
				Q += B(i,j) * s_i * s_j;
			}
		}
		
		Q /= (4.0 * m_TotalNumberOfEdges);
		return Q;
	}
	
} // end namespace anima


