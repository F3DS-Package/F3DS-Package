# Fortran Fine-volume Fluid Dynamics Solver (F3DS) Flamework & Collection

FC=gfortran
#FCFLAGS=-O0 -g -pg -ffree-line-length-none -cpp -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -D_DEBUG # for debug build (gfortran)
FCFLAGS=-O3 -march=native -ffree-line-length-none -fopenmp -cpp # for release build (gfortran)
#FCFLAGS=-g -check all -fpe0 -warn -traceback -debug extended # for debug build (ifort)
#FCFLAGS=-fast -qopenmp

OBJDIR=objs
MODDIR=mods
BINDIR=bins

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
# Tier 3 layer (Utility & Abstract)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/math/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/utils/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/abstracts/*.f90)
# Tier 2 layer (Scheme, cellsystem, and I/O)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/cellsystem/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/eos/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/reconstruction/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/time_stepping/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/grid_parser/*.f90)
FRAMEWORK_SRCS+=$(wildcard $(FRAMEWORK_SRCDIR)/initial_condition_parser/*.f90)
# Tier 1 layer (Application)

# F3DS Collection
# Five-equation model
FIVE_EQ_MODEL_SRCDIR=collection/five_equation_model/src
FIVE_EQ_MODEL_SRCS=$(wildcard $(FIVE_EQ_MODEL_SRCDIR)/nonviscosity_flux/*.f90)
FIVE_EQ_MODEL_SRCS+=$(wildcard $(FIVE_EQ_MODEL_SRCDIR)/space/*.f90)
FIVE_EQ_MODEL_SRCS+=$(wildcard $(FIVE_EQ_MODEL_SRCDIR)/*.f90)

TARGET=mfs.exe
OBJS=$(subst $(JSONFORTDIR)/,$(OBJDIR)/, $(JSONFORTSRCS))
OBJS+=$(subst $(PENF_DIR)/, $(OBJDIR)/, $(PENF_SRC))
OBJS+=$(subst $(FACE_DIR)/, $(OBJDIR)/, $(FACE_SRC))
OBJS+=$(subst $(BEFOR64_DIR)/, $(OBJDIR)/, $(BEFOR64_SRC))
OBJS+=$(subst $(STRINGIFOR_DIR)/, $(OBJDIR)/, $(STRINGIFOR_SRC))
OBJS+=$(subst $(FoXy_DIR)/, $(OBJDIR)/, $(FoXy_SRC))
OBJS+=$(subst $(VTKFORTRAN_DIR)/, $(OBJDIR)/, $(VTKFORTRAN_SRC))
OBJS+=$(subst $(FRAMEWORK_SRCDIR)/, $(OBJDIR)/, $(FRAMEWORK_SRCS))
OBJS+=$(subst $(FIVE_EQ_MODEL_SRCDIR)/, $(OBJDIR)/, $(FIVE_EQ_MODEL_SRCS))
OBJS:=$(subst .f90,.o,$(OBJS))
OBJS:=$(subst .F90,.o,$(OBJS))

$(TARGET): $(OBJS)
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

# Five-equation model
$(OBJDIR)/%.o: $(FIVE_EQ_MODEL_SRCDIR)/%.f90
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	[ -d $(MODDIR) ] || mkdir -p $(MODDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)
	rm -rf $(TARGET)
	rm -rf $(BINDIR)

list:
	echo $(OBJS)
