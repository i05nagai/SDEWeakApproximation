#Suffix Rule
.SUFFIXES :
.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) $(INC_FLAGS) -c $<

include ../common.mk

SDE_WA_SRC =  $(wildcard *.c)
SDE_WA_OBJS =  $(SDE_WA_SRC:.c=.o)
TARGET_SDE_WA = libmt19937.a
TARGETS = $(TARGET_SDE_WA) 

.PHONY: all
all: $(SDE_WA_OBJS)

$(TARGET_SDE_WA): $(SDE_WA_OBJS)
	rm -f $@
	ar cvr $@ $(SDE_WA_OBJS)

Makefile : $(SELF)
	rm -f $@
	cp $(SELF) $@
	chmod +w $@
	echo '# Automatically-generated dependencies list:' >>$@
	gcc ${CFLAGS} ${INC_FLAGS} -MM	\
	${SDE_WA_SRC}	\
	>> $@
	chmod -w $@

.PHONY: clean
clean :
	rm -f *.o $(TARGETS)

