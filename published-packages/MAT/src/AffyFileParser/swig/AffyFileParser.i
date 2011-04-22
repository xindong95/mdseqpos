%module AffyFileParser 
%include "stl.i"

%{
#include "affy-base-types.h"
#include "IntervalEntry.h"
#include "TagValuePair.h"
#include "CDFFileData.h"
#include "CHPFileData.h"
#include "CHPFileWriter.h"
#include "BEDFileData.h"
#include "BEDFileWriter.h"
#include "BARFileData.h"
#include "BARFileWriter.h"
#include "BPMAPFileData.h"
#include "BPMAPFileWriter.h"
#include "MSKFileData.h"
#include "MSKFileWriter.h"
#include "CELFileData.swig.h"
#undef SWIG_croak
#define SWIG_croak(x) { puts(x); SWIG_SetError(x); goto fail; }
%}

%rename (_TagValuePairType_assign) operator=(_TagValuePairType vp);
%rename (_TagValuePairType_equals_obj) operator==(_TagValuePairType vp);
%rename (_TagValuePairType_equals_tag) operator==(const char *tag);

%rename (_BackgroundZoneType_assign) operator=(_BackgroundZoneType zn);
%rename (CExpressionProbeSetResults_assign) operator=(CExpressionProbeSetResults &src);
%rename (CGenotypeProbeSetResults_assign) operator=(CGenotypeProbeSetResults &src);
%rename (CUniversalProbeSetResults_assign) operator=(CUniversalProbeSetResults &src);
%rename (GetQCProbeSetInformation_By_Type) GetQCProbeSetInformation(GeneChipQCProbeSetType qcType, CCDFQCProbeSetInformation & info);

%rename (_GDACSequenceHitItemType_less_than) operator<(const _GDACSequenceHitItemType &rhs) const;
%rename (CGDACSequenceItemWriter_less_than) operator<(const CGDACSequenceItemWriter &rhs) const;

%include "std_string.i"
%apply const std::string& {std::string *};
%include "std_vector.i"

%include "affy-base-types.h"
%include "IntervalEntry.h"
%include "TagValuePair.h"
%include "CDFFileData.h"
%include "CHPFileData.h"
%include "CHPFileWriter.h"
%include "BEDFileData.h"
%include "BEDFileWriter.h"
%include "BARFileData.h"
%include "BARFileWriter.h"
%include "BPMAPFileData.h"
%include "BPMAPFileWriter.h"
%include "MSKFileData.h"
%include "MSKFileWriter.h"
%include "CELFileData.swig.h"
