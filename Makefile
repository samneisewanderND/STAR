# user may define the whole sets of flags
# LDFLAGS
# CPPFLAGS
# CXXFLAGS
# CFLAGS

# or these user-set flags that will be added to standard flags
LDFLAGSextra ?=
CXXFLAGSextra ?=

# user may define the compiler
CXX ?= g++

# project layout
SRCDIR   := src
INCDIR   := include
BUILDDIR := build
BINDIR   := bin
EXTDIR   := external

# vendored third-party libraries
HTSLIB_DIR := $(EXTDIR)/htslib
HTSLIB_LIB := $(HTSLIB_DIR)/libhts.a
OPAL_DIR   := $(EXTDIR)/opal
OPAL_OBJ   := $(BUILDDIR)/opal/opal.o

# include/ holds project headers; src/ is still needed for VERSION;
# external/ provides the htslib/, opal/ and SimpleGoodTuring/ headers;
# build/ holds the generated parametersDefault.xxd
CPPFLAGS_include := -I$(INCDIR) -I$(SRCDIR) -I$(EXTDIR) -I$(BUILDDIR)

# generate a .d dependency file next to each object
DEPFLAGS := -MMD -MP

# pre-defined flags
LDFLAGS_shared := -pthread -L$(HTSLIB_DIR) -Bstatic -lhts -Bdynamic -lz
LDFLAGS_static := -static -static-libgcc -pthread -L$(HTSLIB_DIR) -lhts -lz
LDFLAGS_Mac :=-pthread -lz $(HTSLIB_LIB)
LDFLAGS_Mac_static :=-pthread -lz -static-libgcc $(HTSLIB_LIB)
LDFLAGS_gdb := $(LDFLAGS_shared)

DATE_FMT = --iso-8601=seconds
ifdef SOURCE_DATE_EPOCH
    BUILD_DATE ?= $(shell date -u -d "@$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u -r "$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u "$(DATE_FMT)")
else
    BUILD_DATE ?= $(shell date "$(DATE_FMT)")
endif

BUILD_PLACE ?= $(HOSTNAME):$(shell pwd)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(BUILD_DATE) $(BUILD_PLACE)"'


GIT_CHECK := $(shell git status 1> /dev/null 2> /dev/null && echo 0)
ifeq ($(GIT_CHECK),0)
	GIT_DIFF := $(shell git diff --name-only | tr '\n' ' ')
    GIT_BRANCH_COMMIT_DIFF := $(shell git status | head -n 1) ; \
                                                    $(shell git log | head -n 1) ; \
                                                    diff files: $(GIT_DIFF)
endif
GIT_BRANCH_COMMIT_DIFF := -D'GIT_BRANCH_COMMIT_DIFF="$(GIT_BRANCH_COMMIT_DIFF)"'

# Defaults, can be overridden by make arguments or environment
CXXFLAGS ?= -pipe -Wall -Wextra
CFLAGS ?= -pipe -Wall -Wextra -O3
CXXFLAGS_SIMD ?= -mavx2

# Unconditionally set essential flags and optimization options
CXXFLAGS_common := -std=c++11 -fopenmp $(COMPTIMEPLACE) $(GIT_BRANCH_COMMIT_DIFF)
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS_gdb := -O0 -g3 $(CXXFLAGS_common)

##########################################################################################################
OBJECT_NAMES = systemFunctions.o funPrimaryAlignMark.o \
	SoloFeature_collapseUMI_Graph.o SoloFeature_collapseUMIall.o ParametersClip_initialize.o ClipMate_clip.o ClipCR4.o ClipMate_clipChunk.o ClipMate_initialize.o \
	SoloFeature_loadRawMatrix.o SoloFeature_emptyDrops_CR.o soloInputFeatureUMI.o SoloFeature_countSmartSeq.o SoloFeature_redistributeReadsByCB.o \
	SoloFeature_quantTranscript.o SoloFeature_sumThreads.o SoloFeature_countVelocyto.o SoloFeature_countCBgeneUMI.o \
	Transcriptome_classifyAlign.o Transcriptome_geneFullAlignOverlap_ExonOverIntron.o Transcriptome_alignExonOverlap.o \
	SoloFeature_cellFiltering.o \
	SoloFeature_statsOutput.o bamSortByCoordinate.o SoloBarcode.o \
	ParametersSolo.o SoloRead.o SoloRead_record.o \
	SoloReadBarcode.o SoloReadBarcode_getCBandUMI.o SoloBarcode_extractBarcode.o \
	SoloReadFeature.o SoloReadFeature_record.o SoloReadFeature_inputRecords.o \
	Solo.o SoloFeature.o SoloFeature_outputResults.o SoloFeature_processRecords.o SoloFeature_addBAMtags.o \
	ReadAlign_transformGenome.o Genome_transformGenome.o Transcript_convertGenomeCigar.o \
	twoPassRunPass1.o samHeaders.o Genome_genomeLoad.o Genome_genomeOutLoad.o Transcript_transformGenome.o ReadAlign_outputSpliceGraphSAM.o \
	ReadAlign_mapOneReadSpliceGraph.o SpliceGraph.o SpliceGraph_swScoreSpliced.o SpliceGraph_swTraceBack.o \
	SpliceGraph_findSuperTr.o sjAlignSplit.o \
	GTF.o GTF_transcriptGeneSJ.o GTF_superTranscript.o SuperTranscriptome.o \
	ReadAlign_outputAlignments.o  \
	ReadAlign.o STAR.o \
	SharedMemory.o PackedArray.o SuffixArrayFuns.o Parameters.o Parameters_samAttributes.o InOutStreams.o SequenceFuns.o Genome.o ParametersGenome.o Stats.o \
	Transcript.o Transcript_alignScore.o Transcript_generateCigarP.o Chain.o \
	Transcript_variationAdjust.o Variation.o ReadAlign_waspMap.o \
	ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
	ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o ReadAlign_calcCIGAR.o \
	ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o \
	ReadAlign_peOverlapMergeMap.o ReadAlign_mappedFilter.o \
	ParametersChimeric_initialize.o ReadAlign_chimericDetection.o ReadAlign_chimericDetectionOld.o ReadAlign_chimericDetectionOldOutput.o\
	ChimericDetection.o ChimericDetection_chimericDetectionMult.o ReadAlign_chimericDetectionPEmerged.o \
	stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o \
	ChimericSegment.o ChimericAlign.o ChimericAlign_chimericJunctionOutput.o ChimericAlign_chimericBAMoutput.o ChimericAlign_chimericStitching.o \
	Genome_genomeGenerate.o genomeParametersWrite.o genomeScanFastaFiles.o genomeSAindex.o \
	Genome_insertSequences.o insertSeqSA.o funCompareUintAndSuffixes.o funCompareUintAndSuffixesMemcmp.o \
	TimeFunctions.o ErrorWarning.o streamFuns.o stringSubstituteAll.o \
	Transcriptome.o Transcriptome_quantAlign.o Transcriptome_geneFullAlignOverlap.o \
	ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
	sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o mapThreadsSpawn.o \
	Parameters_readFilesInit.o Parameters_openReadsFiles.o Parameters_closeReadsFiles.o Parameters_readSAMheader.o \
	bam_cat.o serviceFuns.o GlobalVariables.o \
	BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o

OBJECTS := $(addprefix $(BUILDDIR)/,$(OBJECT_NAMES)) $(OPAL_OBJ)

# objects include htslib headers, so build htslib first (order-only)
$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp | $(HTSLIB_LIB)
	@mkdir -p $(@D)
	@echo Compiling $@...
	@$(CXX) -c $(DEPFLAGS) $(CPPFLAGS_include) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(BUILDDIR)/%.o : $(SRCDIR)/%.c | $(HTSLIB_LIB)
	@mkdir -p $(@D)
	$(CXX) -c $(DEPFLAGS) $(CPPFLAGS_include) $(CPPFLAGS) $(CFLAGS) $< -o $@

all: cleanCompileInfo STAR$(SFX)

$(OPAL_OBJ) : $(OPAL_DIR)/opal.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(DEPFLAGS) -I$(OPAL_DIR) -std=c++11 $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGSextra) $(CXXFLAGS_SIMD) $< -o $@

$(BINDIR):
	mkdir -p $@

.PHONY: clean
clean:
	'rm' -rf $(BUILDDIR) $(BINDIR)

.PHONY: CLEAN
CLEAN: clean
	$(MAKE) -C $(HTSLIB_DIR) clean


.PHONY: clean_solo
clean_solo:
	'rm' -f $(BUILDDIR)/Solo*.o $(BUILDDIR)/Solo*.d

.PHONY: cleanCompileInfo
cleanCompileInfo:
	'rm' -f $(BUILDDIR)/STAR.o $(BUILDDIR)/Parameters.o

-include $(OBJECTS:.o=.d)

.PHONY: htslib
htslib : $(HTSLIB_LIB)

$(HTSLIB_LIB) :
	$(MAKE) -C $(HTSLIB_DIR) lib-static

# run xxd from inside source/ so the generated array is named parametersDefault
$(BUILDDIR)/parametersDefault.xxd: $(SRCDIR)/parametersDefault
	@mkdir -p $(@D)
	cd $(SRCDIR) && xxd -i parametersDefault > $(CURDIR)/$@

$(BUILDDIR)/STAR.o $(BUILDDIR)/Parameters.o: $(BUILDDIR)/parametersDefault.xxd

# link targets are named after the build variant, not the output file
.PHONY: STAR$(SFX) STARstatic$(SFX) STARlong$(SFX) STARlongStatic$(SFX) POSIXSHARED gdb gdb-long STARforMacStatic STARlongForMacStatic

STAR$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STAR$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STAR$(SFX) : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARstatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STARstatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARstatic$(SFX) : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlong$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlong$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STARlong$(SFX) : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlongStatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARlongStatic$(SFX) : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)



POSIXSHARED : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -DPOSIX_SHARED_MEM $(CXXFLAGS)
POSIXSHARED : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
POSIXSHARED : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) $(CXXFLAGS)
gdb : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS)
gdb : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb-long : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
gdb-long : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS)
gdb-long : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARforMacStatic : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARforMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARforMacStatic : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongForMacStatic : CXXFLAGS := -D'COMPILE_FOR_LONG_READS' $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARlongForMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARlongForMacStatic : $(OBJECTS) $(HTSLIB_LIB) | $(BINDIR)
	$(CXX) -o $(BINDIR)/STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)
