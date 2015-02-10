#
# makefile for libsde_wa.a
#
# $Source$
# $Author$
# $Date$
# $Revision$
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; version 2.1 of the License.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License version 2.1 for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# version 2.1 along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, US

#Suffix Rule
.SUFFIXES :
.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) $(INC_FLAGS) -c $<

include ../common.mk

SDE_WA_SRC =  $(wildcard *.c)
SDE_WA_OBJS =  $(SDE_WA_SRC:.c=.o)
TARGET_SDE_WA = libsde_wa.a
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

