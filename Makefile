# Fortran Finite-volume Fluid Dynamics Solver (F3DS) Flamework & Collection

# User inputs
COMPILER = gfortran
DEBUG    = no

# Constants
GNU_RELEAS_FLAGS =-O3 -march=native -ffree-line-length-none -cpp -fopenmp
GNU_DEBUG_FLAGS  =-O0 -g -pg -ffree-line-length-none -cpp -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=nan -D_DEBUG

INTEL_RELEAS_FLAGS =-fast -fpp -std03 -O3 -ipo -inline all -ipo-jobs4 -vec-report1 -qopenmp
INTEL_DEBUG_FLAGS  =-g -std03  -check all -fpe0 -warn -traceback -debug extended -fpp -D_DEBUG

OBJDIR=objs
MODDIR=mods
BINDIR=bins

# Setup
FC=$(COMPILER)
ifeq "$(DEBUG)" "no"
	FCFLAGS=$(if $(findstring gfortran, $(FC)), $(GNU_RELEAS_FLAGS), $(INTEL_RELEAS_FLAGS))
else
	FCFLAGS=$(if $(findstring gfortran, $(FC)), $(GNU_DEBUG_FLAGS), $(INTEL_DEBUG_FLAGS))
endif

# Flags for module strage
FCFLAGS+=$(if $(findstring gfortran, $(FC)),-J $(MODDIR),-module $(MODDIR))

# JSON Fortran
JSONFORTDIR=third_party/json-fortran/src
JSONFILES=json_kinds.F90 json_parameters.F90 json_string_utilities.F90 json_value_module.F90 json_file_module.F90 json_module.F90
JSONFORTSRCS=$(foreach file,$(JSONFILES),$(JSONFORTDIR)/$(file))
JSONFORT_OBJS=$(subst $(JSONFORTDIR)/,$(OBJDIR)/, $(JSONFORTSRCS))

# PENF
PENF_DIR=third_party/PENF/src/lib
PENF_FILES=penf_global_parameters_variables.F90 penf_b_size.F90 penf_stringify.F90 penf.F90
PENF_SRC=$(foreach file,$(PENF_FILES),$(PENF_DIR)/$(file))

# FACE
FACE_DIR=third_party/FACE/src/lib
FACE_SRC=$(wildcard $(FACE_DIR)/*.F90)

# BeFoR64
BEFOR64_DIR=third_party/BeFoR64/src/lib
BEFOR64_SRC=$(wildcard $(BEFOR64_DIR)/*.F90)

# StringiFor
STRINGIFOR_DIR=third_party/StringiFor/src/lib
STRINGIFOR_FILES=stringifor_string_t.F90 stringifor.F90
STRINGIFOR_SRC=$(foreach file,$(STRINGIFOR_FILES),$(STRINGIFOR_DIR)/$(file))

# FoXy
FoXy_DIR=third_party/FoXy/src/lib
FoXy_FILES=foxy_xml_tag.F90 foxy_xml_file.f90 foxy.f90
FoXy_SRC=$(foreach file,$(FoXy_FILES),$(FoXy_DIR)/$(file))

# VTKFortran
VTKFORTRAN_DIR=third_party/VTKFortran/src/lib
VTKFORTRAN_FILES=vtk_fortran_parameters.f90 vtk_fortran_dataarray_encoder.f90 \
vtk_fortran_vtk_file_xml_writer_abstract.f90 \
vtk_fortran_vtk_file_xml_writer_appended.f90 \
vtk_fortran_vtk_file_xml_writer_ascii_local.f90 \
vtk_fortran_vtk_file_xml_writer_binary_local.f90 \
vtk_fortran_pvtk_file.f90 \
vtk_fortran_vtk_file.f90 vtk_fortran_vtm_file.F90 vtk_fortran.f90
VTKFORTRAN_SRC=$(foreach file,$(VTKFORTRAN_FILES),$(VTKFORTRAN_DIR)/$(file))

# F3DS Flamework
FRAMEWORK_SRCDIR=framework/src
# Tier 4 layer (Standart)
FRAMEWORK_SRCS=$(wildcard  $(FRAMEWORK_SRCDIR)/std/*.f90)
# Tier 3 layer (Utility & abstract)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/utils/system_utils/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/utils/math/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/utils/cellsystem_utils/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/utils/reconstructor_utils/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/abstract/auxiliary/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/abstract/core/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/abstract/user_inherited/*.f90)
# Tier 2 layer ()
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/configuration/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/grid_parser/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/initial_condition_parser/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/result_writer/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/eos/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/riemann_solver/hll/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/gradient_calculator/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/interpolator/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/reconstructor/muscl/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/reconstructor/weno/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/time_stepping/rk/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/time_stepping/explicit/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/termination_criterion/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/time_increment_controller/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/parallelizer/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/measurement_tools/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/face_gradient_interpolator/*.f90)
# Tier 1 layer
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/generator/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/cellsystem/*.f90)

# Flamework objects
FRAMEWORK_OBJS=$(subst $(JSONFORTDIR)/,$(OBJDIR)/, $(JSONFORTSRCS))
FRAMEWORK_OBJS+=$(subst $(PENF_DIR)/, $(OBJDIR)/, $(PENF_SRC))
FRAMEWORK_OBJS+=$(subst $(FACE_DIR)/, $(OBJDIR)/, $(FACE_SRC))
FRAMEWORK_OBJS+=$(subst $(BEFOR64_DIR)/, $(OBJDIR)/, $(BEFOR64_SRC))
FRAMEWORK_OBJS+=$(subst $(STRINGIFOR_DIR)/, $(OBJDIR)/, $(STRINGIFOR_SRC))
FRAMEWORK_OBJS+=$(subst $(FoXy_DIR)/, $(OBJDIR)/, $(FoXy_SRC))
FRAMEWORK_OBJS+=$(subst $(VTKFORTRAN_DIR)/, $(OBJDIR)/, $(VTKFORTRAN_SRC))
FRAMEWORK_OBJS+=$(subst $(FRAMEWORK_SRCDIR)/, $(OBJDIR)/, $(FRAMEWORK_SRCS))
FRAMEWORK_OBJS:=$(subst .f90,.o,$(FRAMEWORK_OBJS))
FRAMEWORK_OBJS:=$(subst .F90,.o,$(FRAMEWORK_OBJS))

# F3DS Collection
# Five-equation model common tools
FIVE_EQUATION_MODEL_COMMON_SRCDIR=collection/five_equation_model_common/src
FIVE_EQUATION_MODEL_COMMON_SRCS=$(wildcard  $(FIVE_EQUATION_MODEL_COMMON_SRCDIR)/special_reconstructor/*.f90)
FIVE_EQUATION_MODEL_COMMON_SRCS+=$(wildcard $(FIVE_EQUATION_MODEL_COMMON_SRCDIR)/special_generator/*.f90)
FIVE_EQUATION_MODEL_COMMON_OBJS=$(FRAMEWORK_OBJS)
FIVE_EQUATION_MODEL_COMMON_OBJS+=$(subst .f90,.o, $(subst $(FIVE_EQUATION_MODEL_COMMON_SRCDIR)/, $(OBJDIR)/, $(FIVE_EQUATION_MODEL_COMMON_SRCS)))

# Five-equation model
FIVE_EQUATION_MODEL_TARGET=f5eq
FIVE_EQUATION_MODEL_SRCDIR=collection/five_equation_model/src
FIVE_EQUATION_MODEL_SRCS=$(wildcard $(FIVE_EQUATION_MODEL_SRCDIR)/variables/*.f90)
FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(FIVE_EQUATION_MODEL_SRCDIR)/model/*.f90)
FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(FIVE_EQUATION_MODEL_SRCDIR)/special_reconstructor/*.f90)
FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(FIVE_EQUATION_MODEL_SRCDIR)/special_generator/*.f90)
FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(FIVE_EQUATION_MODEL_SRCDIR)/*.f90)
FIVE_EQUATION_MODEL_OBJS=$(FIVE_EQUATION_MODEL_COMMON_OBJS)
FIVE_EQUATION_MODEL_OBJS+=$(subst .f90,.o, $(subst $(FIVE_EQUATION_MODEL_SRCDIR)/, $(OBJDIR)/, $(FIVE_EQUATION_MODEL_SRCS)))

# Viscous five-equation model
VISCOUS_FIVE_EQUATION_MODEL_TARGET=fv5eq
VISCOUS_FIVE_EQUATION_MODEL_SRCDIR=collection/viscous_five_equation_model/src
VISCOUS_FIVE_EQUATION_MODEL_SRCS=$(wildcard $(VISCOUS_FIVE_EQUATION_MODEL_SRCDIR)/variables/*.f90)
VISCOUS_FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(VISCOUS_FIVE_EQUATION_MODEL_SRCDIR)/model/*.f90)
VISCOUS_FIVE_EQUATION_MODEL_SRCS+=$(wildcard $(VISCOUS_FIVE_EQUATION_MODEL_SRCDIR)/*.f90)
VISCOUS_FIVE_EQUATION_MODEL_OBJS=$(FIVE_EQUATION_MODEL_COMMON_OBJS)
VISCOUS_FIVE_EQUATION_MODEL_OBJS+=$(subst .f90,.o, $(subst $(VISCOUS_FIVE_EQUATION_MODEL_SRCDIR)/, $(OBJDIR)/, $(VISCOUS_FIVE_EQUATION_MODEL_SRCS)))

all: $(FIVE_EQUATION_MODEL_TARGET) $(VISCOUS_FIVE_EQUATION_MODEL_TARGET)

# Help
help:
	@echo -e '\033[1;32m Make options of F3DS framework & collection \033[0m'
	@echo -e '\033[1;32m Usage: \033[0m'
	@echo -e '  make {options}'
	@echo -e '\033[1;32m Options: \033[0m'
	@echo -e '  COMPILER={compiler name} : Compiler name you want to use. F3DS support for ifort and gfortran now!'
	@echo -e '  DEBUG={yes/no}           : If you set DEBUG=yes, Debug infomations are embeded to your code.'
	@echo -e '\033[1;32m Example use: \033[0m'
	@echo -e '  make COMPILER=gfortran DEBUG=no'

# Five-equation model compiling rule
$(FIVE_EQUATION_MODEL_TARGET): $(FIVE_EQUATION_MODEL_OBJS)
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) $+ -o $(BINDIR)/$@

# Viscous five-equation model compiling rule
$(VISCOUS_FIVE_EQUATION_MODEL_TARGET): $(VISCOUS_FIVE_EQUATION_MODEL_OBJS)
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) $+ -o $(BINDIR)/$@

# JSON Fortran
$(OBJDIR)/%.o: $(JSONFORTDIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# PENF
$(OBJDIR)/%.o: $(PENF_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# FACE
$(OBJDIR)/%.o: $(FACE_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# BeFoR64
$(OBJDIR)/%.o: $(BEFOR64_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# StringiFor
$(OBJDIR)/%.o: $(STRINGIFOR_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# FoXy
$(OBJDIR)/%.o: $(FoXy_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@
$(OBJDIR)/%.o: $(FoXy_DIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# VTKFortran
$(OBJDIR)/%.o: $(VTKFORTRAN_DIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@
$(OBJDIR)/%.o: $(VTKFORTRAN_DIR)/%.F90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# F3DS Flamework
$(OBJDIR)/%.o: $(FRAMEWORK_SRCDIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Five-equation model common tools
$(OBJDIR)/%.o: $(FIVE_EQUATION_MODEL_COMMON_SRCDIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Five-equation model
$(OBJDIR)/%.o: $(FIVE_EQUATION_MODEL_SRCDIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Viscous five-equation model
$(OBJDIR)/%.o: $(VISCOUS_FIVE_EQUATION_MODEL_SRCDIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)
	rm -rf $(BINDIR)

list:
	@echo -e '\033[1;32m Framework objects:\033[0m'
	@echo -e $(FRAMEWORK_OBJS)
	@echo -e '\033[1;32m Viscous 5-equation model objects:\033[0m'
	@echo -e $(VISCOUS_FIVE_EQUATION_MODEL_OBJS)
	@echo -e '\033[1;32m 5-equation model objects objects:\033[0m'
	@echo -e $(FIVE_EQUATION_MODEL_OBJS)
