// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIgiacomodIcaendgzmIsipmanalysisdImacrodI04_multi_scandIsipm_tot_analysis_cpp_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/giacomo/caendgz-sipmanalysis/macro/04_multi_scan/./sipm_tot_analysis.cpp"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_sipm_tot_analysis_cpp_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./sipm_tot_analysis.cpp",
nullptr
    };
    static const char* includePaths[] = {
"/home/giacomo/tools/root/include",
"/home/giacomo/tools/root/etc/",
"/home/giacomo/tools/root/etc//cling",
"/home/giacomo/tools/root/etc//cling/plugins/include",
"/home/giacomo/tools/root/include/",
"/home/giacomo/tools/root/include",
"/home/giacomo/caendgz-sipmanalysis/macro/04_multi_scan/",
"/home/giacomo/tools/root/include/",
"/home/giacomo/caendgz-sipmanalysis/macro/04_multi_scan/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "sipm_tot_analysis_cpp_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "sipm_tot_analysis_cpp_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./sipm_tot_analysis.cpp"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"findCalibratedCutoffsInteractive", payloadCode, "@",
"g_analysis_mode", payloadCode, "@",
"sipm_tot_analysis", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("sipm_tot_analysis_cpp_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_sipm_tot_analysis_cpp_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_sipm_tot_analysis_cpp_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_sipm_tot_analysis_cpp_ACLiC_dict() {
  TriggerDictionaryInitialization_sipm_tot_analysis_cpp_ACLiC_dict_Impl();
}
