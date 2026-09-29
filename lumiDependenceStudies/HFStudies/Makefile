CXX = g++
# CXXFLAGS += -O2 -Wall -Werror -Wextra
ROOTFLAGS := `root-config --cflags --libs` -I../../include/

BUILDDIR = ./build

SRCS  = $(wildcard *.cc)
EXES  = $(patsubst %.cc,%.exe,$(SRCS))
DEPS  = $(patsubst %.cc,$(BUILDDIR)/%.d,$(SRCS))

.PHONY: all clean

all: $(EXES)

%.exe: %.cc
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(ROOTFLAGS) $(CXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

clean:
	@$(RM) $(EXES) $(DEPS)
	@[ -n "$(BUILDDIR)" ] && rm -rf "$(BUILDDIR)"

-include $(DEPS)
