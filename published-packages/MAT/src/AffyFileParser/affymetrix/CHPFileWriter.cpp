/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

#include "CHPFileWriter.h"
#include "FileWriter.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <istream>
#include <fstream>
#include <cstring>		// port to gcc 4.3, by Tao Liu
#include <climits>		// port to gcc 4.3, by Tao Liu
#include <cstdlib>		// port to gcc 4.3, by Tao Liu
#include <cstdio>		// port to gcc 4.3, by Tao Liu
#include <cctype>		// port to gcc 4.3, by Tao Liu
#include <algorithm>		// port to gcc 4.3, by Tao Liu
#include <iterator>		// port to gcc 4.3, by Tao Liu
#include <memory>		// port to gcc 4.3, by Tao Liu
#include <typeinfo>		// port to gcc 4.3, by Tao Liu

#ifdef _INCLUDE_UNISTD_HEADER_
#include <unistd.h>
#endif

using namespace affxchp;
using namespace affxchpwriter;

//////////////////////////////////////////////////////////////////////

#define DELIMCHAR 0x14
#define MIN_CELLSTR 4
#define CHP_FILE_MAGIC_NUMBER 65
#define CHP_FILE_VERSION_NUMBER 2
#define EXPRESSION_ABSOLUTE_STAT_ANALYSIS 2
#define EXPRESSION_COMPARISON_STAT_ANALYSIS 3

//////////////////////////////////////////////////////////////////////

CCHPFileWriter::CCHPFileWriter() : CCHPFileData()
{
}

//////////////////////////////////////////////////////////////////////

CCHPFileWriter::~CCHPFileWriter()
{
	Clear();
}

//////////////////////////////////////////////////////////////////////

bool CCHPFileWriter::CreateNewFile()
{
	// Open the file.
	m_strError = "";
	m_NewChpFile.open(m_FileName.c_str(), std::ios::out | std::ios::binary);
	if (!m_NewChpFile)
	{
		m_strError = "Unable to open the file.";
		return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::InitializeForWriting(affxcdf::CCDFFileData &cdf)
{
	InitializeForWriting(
		cdf.GetHeader().GetRows(),
		cdf.GetHeader().GetCols(),
		cdf.GetHeader().GetNumProbeSets(),
		cdf.GetChipType().c_str(),
		cdf.GetProbeSetType(0));
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::InitializeForWriting(int numRows, int numCols, int numProbeSets, const char *chipType, affxcdf::GeneChipProbeSetType probeSetType)
{
	// Set the header values.
	m_Header.SetCols(numCols);
	m_Header.SetRows(numRows);
	m_Header.SetNumProbeSets(numProbeSets);
	m_Header.SetChipType(chipType);

	switch (probeSetType)
	{
	case affxcdf::UnknownProbeSetType:
		m_Header.SetAssayType(CCHPFileHeader::Unknown);
		break;

	case affxcdf::ExpressionProbeSetType:
		m_Header.SetAssayType(CCHPFileHeader::Expression);
		break;

	case affxcdf::GenotypingProbeSetType:
		m_Header.SetAssayType(CCHPFileHeader::Genotyping);
		break;

	case affxcdf::ResequencingProbeSetType:
		m_Header.SetAssayType(CCHPFileHeader::Resequencing);
		break;

	case affxcdf::TagProbeSetType:
		m_Header.SetAssayType(CCHPFileHeader::Universal);
		break;

	default:
		m_Header.SetAssayType(CCHPFileHeader::Unknown);
		break;
	}


	// Allocate memory for probe set results
	m_ProbeSetResults.resize(numProbeSets);
	CProbeSetResults *pResults;
	for (int iset=0; iset<numProbeSets; iset++)
	{
		switch (probeSetType)
		{
			case affxcdf::ExpressionProbeSetType:
				pResults = new CExpressionProbeSetResults;
				break;

			case affxcdf::GenotypingProbeSetType:
				pResults = new CGenotypeProbeSetResults;
				break;

			default:
				pResults = NULL;
				break;
		}
		m_ProbeSetResults[iset] = pResults;
	}
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetParentCelFileName(const char *str)
{
	m_Header.SetParentCellFile(str);
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetProgID(const char *str)
{
	m_Header.SetProgID(str);
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetAlgorithmName(const char *str)
{
	std::string name = str;
	if (m_Header.GetAssayType() == CCHPFileHeader::Expression)
	{
		name = "Expression" + std::string(str);
	}

	else if (m_Header.GetAssayType() == CCHPFileHeader::Genotyping)
	{
		name = "Genotyping" + std::string(str);
	}

	m_Header.SetAlgName(name.c_str());
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetAlgorithmVersion(const char *str)
{
	m_Header.SetAlgVersion(str);
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::AddAlgorithmParameter(const char *tag, const char *value)
{
	TagValuePairType param;
	param.Tag = tag;
	param.Value = value;
	m_Header.AlgorithmParameters().push_back(param);
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::AddChipSummaryParameter(const char *tag, const char *value)
{
	TagValuePairType param;
	param.Tag = tag;
	param.Value = value;
	m_Header.SummaryParameters().push_back(param);
}


//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::AddBackgroundInfo(int nZones, float smoothFactor)
{
	m_Header.GetBackgroundZoneInfo().number_zones = nZones;
	m_Header.GetBackgroundZoneInfo().smooth_factor = smoothFactor;
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::AddBackgroundZone(int x, int y, float value)
{
	BackgroundZoneType zone;
	zone.centerx = (float) x;
	zone.centery = (float) y;
	zone.background = value;
	m_Header.GetBackgroundZoneInfo().zones.push_back(zone);
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetExpressionEntry(int index, CExpressionProbeSetResults *pEntry)
{
	CExpressionProbeSetResults *pExpressionResult = GetExpressionResults(index);
	*pExpressionResult = *pEntry;
}

//////////////////////////////////////////////////////////////////////

void CCHPFileWriter::SetMappingEntry(int index, CGenotypeProbeSetResults *pEntry)
{
	CGenotypeProbeSetResults *pMappingResult = GetGenotypingResults(index);
	*pMappingResult = *pEntry;
}

//////////////////////////////////////////////////////////////////////

bool CCHPFileWriter::Save()
{
	// Only continue if genotyping or expression
	if (m_Header.GetAssayType() != CCHPFileHeader::Expression && 
		m_Header.GetAssayType() != CCHPFileHeader::Genotyping)
	{
		m_strError = "The software only supports writing expression or genotyping CHP files.";
		return false;
	}

	
	// Write the header
	int magic = CHP_FILE_MAGIC_NUMBER;
	WRITE_INT(m_NewChpFile, magic);
	int version = 1;
	WRITE_INT(m_NewChpFile, version);
	WRITE_USHORT(m_NewChpFile, m_Header.GetCols());
	WRITE_USHORT(m_NewChpFile, m_Header.GetRows());
	WRITE_INT(m_NewChpFile, m_Header.GetNumProbeSets());
	WRITE_INT(m_NewChpFile, 0); // no qc data extracted.
	WRITE_INT(m_NewChpFile, m_Header.GetAssayType());
	WRITE_STRING(m_NewChpFile, m_Header.GetProgID());
	WRITE_STRING(m_NewChpFile, m_Header.GetParentCellFile());
	WRITE_STRING(m_NewChpFile, m_Header.GetChipType());
	WRITE_STRING(m_NewChpFile, m_Header.GetAlgName());
	WRITE_STRING(m_NewChpFile, m_Header.GetAlgVersion());
	WRITE_INT(m_NewChpFile, (int) m_Header.AlgorithmParameters().size());
	TagValuePairTypeList::iterator iter;
	for (iter=m_Header.AlgorithmParameters().begin(); iter!=m_Header.AlgorithmParameters().end(); ++iter)
	{
		WRITE_STRING(m_NewChpFile, iter->Tag);
		WRITE_STRING(m_NewChpFile, iter->Value);
	}
	WRITE_INT(m_NewChpFile, (int) m_Header.SummaryParameters().size());
	for (iter=m_Header.SummaryParameters().begin(); iter!=m_Header.SummaryParameters().end(); ++iter)
	{
		WRITE_STRING(m_NewChpFile, iter->Tag);
		WRITE_STRING(m_NewChpFile, iter->Value);
	}

	// Write the zone info.
	WRITE_INT(m_NewChpFile, m_Header.GetBackgroundZoneInfo().number_zones);
	WRITE_FLOAT(m_NewChpFile, m_Header.GetBackgroundZoneInfo().smooth_factor);
	BackgroundZoneTypeList::iterator start(m_Header.GetBackgroundZoneInfo().zones.begin());
	BackgroundZoneTypeList::iterator end(m_Header.GetBackgroundZoneInfo().zones.end());
	BackgroundZoneType zone;
	for (; start != end; ++start)
	{
		zone = (*start);
		WRITE_FLOAT(m_NewChpFile, start->centerx);
		WRITE_FLOAT(m_NewChpFile, start->centery);
		WRITE_FLOAT(m_NewChpFile, start->background);
	}

	// Write the probe set data
	if (m_Header.GetAssayType() == CCHPFileHeader::Expression)
	{
		// Set the type of analysis
		CExpressionProbeSetResults * pResults = GetExpressionResults(0);
		int resultsSize = UCHAR_SIZE + FLOAT_SIZE + FLOAT_SIZE + USHORT_SIZE + USHORT_SIZE; 
		unsigned char analysisType = EXPRESSION_ABSOLUTE_STAT_ANALYSIS;
		if (pResults->m_HasCompResults)
		{
			analysisType = EXPRESSION_COMPARISON_STAT_ANALYSIS;
			resultsSize += UCHAR_SIZE + FLOAT_SIZE + FLOAT_SIZE + FLOAT_SIZE + FLOAT_SIZE + USHORT_SIZE;
		}
		WRITE_UCHAR(m_NewChpFile, analysisType);
		WRITE_INT(m_NewChpFile, resultsSize);

		// Write each probe set result.
		for (int iset=0; iset<m_Header.GetNumProbeSets(); iset++)
		{
			pResults = GetExpressionResults(iset);

			// Write the absolute data.
			WRITE_UCHAR(m_NewChpFile, pResults->Detection);
			WRITE_FLOAT(m_NewChpFile, pResults->DetectionPValue);
			WRITE_FLOAT(m_NewChpFile, pResults->Signal);
			WRITE_USHORT(m_NewChpFile, pResults->NumPairs);
			WRITE_USHORT(m_NewChpFile, pResults->NumUsedPairs);

			// Write the comparison data
			if (pResults->m_HasCompResults == true)
			{
				WRITE_UCHAR(m_NewChpFile, pResults->Change);
				WRITE_FLOAT(m_NewChpFile, pResults->ChangePValue);
				WRITE_FLOAT(m_NewChpFile, pResults->SignalLogRatio);
				WRITE_FLOAT(m_NewChpFile, pResults->SignalLogRatioLow);
				WRITE_FLOAT(m_NewChpFile, pResults->SignalLogRatioHigh);
				WRITE_USHORT(m_NewChpFile, pResults->NumCommonPairs);
			}
		}
	}
	else if (m_Header.GetAssayType() == CCHPFileHeader::Genotyping)
	{
		float fval;
		CGenotypeProbeSetResults * pResults;
		int resultsSize = UCHAR_SIZE + FLOAT_SIZE + FLOAT_SIZE + FLOAT_SIZE + FLOAT_SIZE + FLOAT_SIZE;
		WRITE_INT(m_NewChpFile, resultsSize);
		for (int iset=0; iset<m_Header.GetNumProbeSets(); iset++)
		{
			pResults = GetGenotypingResults(iset);

			// Write probe set result.
			WRITE_UCHAR(m_NewChpFile, pResults->AlleleCall);
			WRITE_FLOAT(m_NewChpFile, pResults->Confidence);

			fval = pResults->pvalue_AA;
			if (fval == 0)
				fval = pResults->RAS1;
			WRITE_FLOAT(m_NewChpFile, fval);

			fval = pResults->pvalue_AB;
			if (fval == 0)
				fval = pResults->RAS2;
			WRITE_FLOAT(m_NewChpFile, fval);

			WRITE_FLOAT(m_NewChpFile, pResults->pvalue_BB);
			WRITE_FLOAT(m_NewChpFile, pResults->pvalue_NoCall);
		}
	}

	// Close the file and check the status.
	m_NewChpFile.close();

	return !m_NewChpFile.fail();
}

//////////////////////////////////////////////////////////////////////
