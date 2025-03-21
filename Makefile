# With help of article: https://habr.com/ru/post/155201/

# Определение операционной системы
ifeq ($(OS),Windows_NT)
	MKDIR = if not exist $(subst /,\,$(1)) mkdir $(subst /,\,$(1))
	RM = if exist $(subst /,\,$(1)) rmdir /s /q $(subst /,\,$(1))
	EXECUTABLE = bin/app.exe
else
	MKDIR = mkdir -p $(1)
	RM = rm -rf $(1)
	EXECUTABLE = bin/app
endif

CC = g++
CFLAGS = -c -Wall -I./include
LDFLAGS =
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	@$(call MKDIR,bin)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

build/%.o: src/%.cpp
	@$(call MKDIR,build)
	$(CC) $(CFLAGS) $< -o $@

clean:
	@$(call RM,build)
	@$(call RM,bin)

.PHONY: all clean
