CC      =	gcc
CFLAGS  =	-Wall -O2  -Wno-unused-function -Wno-misleading-indentation
DFLAGS  =	-g -Wall  
HTSLIB_DIR = ./htslib
HTSLIB  =   $(HTSLIB_DIR)/libhts.a
LIB     =	$(HTSLIB) -lm -lz -lpthread -lbz2 -llzma -lcurl
INCLUDE = -I ./htslib

ifneq ($(gdb),)
	CFLAGS  =	-g -Wall -Wno-unused-variable -Wno-misleading-indentation -Wno-unused-but-set-variable -Wno-unused-function
endif

BIN_DIR =	.
SRC_DIR =   .

HTS_ALL =   hts_all
SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/fxtools

GDB_DEBUG   =   $(BIN_DIR)/gdb_fxtools
NOR_DEBUG   =   $(BIN_DIR)/debug_fxtools
DMARCRO 	=	-D __DEBUG__

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

#all:       $(SOURCE) $(BIN) 
all:       $(HTS_ALL) $(BIN) 
gdb_fxtools: $(SOURCE) $(GDB_DEBUG) 
debug_fxtools: $(SOURCE) $(NOR_DEBUG)


$(HTS_ALL):
	cd $(HTSLIB_DIR); make; cd ../
$(BIN): $(OBJS)
		$(CC) $(OBJS) -o $@ $(LIB)

$(GDB_DEBUG):
		$(CC) $(DFLAGS) $(INCLUDE) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(NOR_DEBUG):
		$(CC) $(CFLAGS) $(INCLUDE) $(SOURCE) $(DMARCRO) -o $@ $(LIB)


clean:
		rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
		rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(NOR_DEBUG)

fxtools.o: fxtools.c fxtools.h
