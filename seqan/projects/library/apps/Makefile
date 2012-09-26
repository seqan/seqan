SEQAN_BASE = ..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif


tbb_root?=../extra/tbb

#check, if tbb_root is not absolute path (the filter keeps only /* paths)
ifeq ($(filter /% $(SLASH)%, $(subst :, ,$(tbb_root)) ),)
    # also changes related variables like work_dir
    override tbb_root := $(CWD)$(SLASH)..
    export TBB21_INSTALL_DIR := $(tbb_root)
endif

# explicitly compile for a 32 or 64 bit platform
#CXXFLAGS += -m32
#CXXFLAGS += -m64

CXXFLAGS += -I$(SEQAN_BASE)
#CXXFLAGS += -I~/Documents/Development/
CXXFLAGS += -O9 -fopenmp -march=native -DSEQAN_ENABLE_DEBUG=0
#CXXFLAGS += -O0 -g -fopenmp -DSEQAN_ENABLE_DEBUG=1
CXXFLAGS += -pedantic -W -Wall -Wno-variadic-macros
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += -DSEQAN_ENABLE_TESTING=0
CXXFLAGS += -static
CXXFLAGS += ${CXXEXTRAFLAGS}

# Intel Threading Building Blocks
#CXXFLAGS += -I$(SEQAN_BASE)/extra/tbb/include/
#LDFLAGS += -L$(SEQAN_BASE)/extra/tbb/build/macos_em64t_gcc_cc4.0.1_os10.5.5_debug/ -ltbb_debug -lpthread

TARGETS = dfi/dfi seqan_tcoffee/seqan_tcoffee seqcons/seqcons razers/paramChooser razers/razers razers2/razers pair_align/pair_align micro_razers/micro_razers tree_recon/tree_recon #transquant/transquant

all:           check_seqan_base $(TARGETS)
dfi:           check_seqan_base dfi/dfi
fiona:         check_seqan_base fiona/fiona
razers:        check_seqan_base razers/razers razers/paramChooser
razers2:       check_seqan_base razers2/razers
razers3:       check_seqan_base razers3/razers
micro_razers:  check_seqan_base micro_razers/micro_razers
seqan_tcoffee: check_seqan_base seqan_tcoffee/seqan_tcoffee
seqcons:       check_seqan_base seqcons/seqcons
pair_align:    check_seqan_base pair_align/pair_align
tree_recon:    check_seqan_base tree_recon/tree_recon
transquant:    check_seqan_base transquant/transquant transquant/rectifier transquant/transsplice
insegt:        check_seqan_base insegt/insegt

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
