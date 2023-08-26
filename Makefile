CC = g++
SRC = src
IDIR = include
OBJDIR = objdir
CFLAGS = -std=gnu++11 -I $(IDIR)
LINK_LIB = -lfftw3f -lfftw3f_omp -fopenmp 

#optimization
O= -O2

VPATH= $(SRC):$(IDIR)

_OBJS= $(patsubst $(SRC)/%.C, %.o, $(wildcard $(SRC)/*.C))
OBJS= $(patsubst %, $(OBJDIR)/%, $(_OBJS))

PowerSpectrum: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LINK_LIB)

$(OBJDIR)/ps_func_mod1.o: ps_func_mod1.C pwr_spc.H
	$(CC) -c -o $@ $< $(CFLAGS) $(O)

$(OBJDIR)/ps_aux.o: ps_aux.C ps_aux.H
	$(CC) -c -o $@ $< $(CFLAGS) $(O)

$(OBJDIR)/ps_main.o: ps_main.C\
	ps_aux.H pwr_spc.H
	$(CC) -c -o $@ $< $(CFLAGS) $(O)

$(OBJDIR)/ps_class.o: ps_class.C pwr_spc_decl.H
	$(CC) -c -o $@ $< $(CFLAGS) $(O)

$(OBJS): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY: clean
clean:
	-rm -r PowerSpectrum $(OBJDIR)
