include ./Makefile.common

.PHONY: sparkbwa libbwa.so bwa clean

all: libbwa.so libcushaw.so
	@echo "================================================================================"
	@echo "SparkBWA library has been built."
	@echo "Location    = $(LOCATION)/$(BUILD_DIR)/"
	@echo "JAVA_HOME   = $(JAVA_HOME)"
	@echo "================================================================================"

bwa:
	$(MAKE) -C $(BWA_DIR)
	if [ ! -d "$(BUILD_DIR)" ]; then mkdir $(BUILD_DIR); fi
	cp $(BWA_DIR)/*.o $(BUILD_DIR)

cushaw:
	$(MAKE) -C $(CUSHAW_DIR)
	if [ ! -d "$(BUILD_DIR)" ]; then mkdir $(BUILD_DIR); fi
	cp $(CUSHAW_DIR)/objs/*.o $(BUILD_DIR)

sparkbwa:
	if [ ! -d "$(BUILD_DIR)" ]; then mkdir $(BUILD_DIR); fi
	$(CC) $(SPARKBWA_FLAGS) bwa_jni.c -o $(BUILD_DIR)/bwa_jni.o $(LIBBWA_LIBS)

sparkcushaw:
	if [ ! -d "$(BUILD_DIR)" ]; then mkdir $(BUILD_DIR); fi
	$(CXX) $(SPARKCUSHAW_FLAGS) cushaw_jni.cpp -o $(BUILD_DIR)/cushaw_jni.o $(LIBCUSHAW_LIBS)

libbwa.so: sparkbwa bwa
	$(CC) $(LIBBWA_FLAGS) $(BUILD_DIR)/libbwa.so $(BUILD_DIR)/*.o $(LIBBWA_LIBS)
	cp $(BUILD_DIR)/libbwa.so ../../../target/classes/
	rm $(BUILD_DIR)/*.o

libcushaw.so: cushaw sparkcushaw
	$(CXX) $(LIBCUSHAW_FLAGS) $(BUILD_DIR)/libcushaw.so $(BUILD_DIR)/*.o $(LIBCUSHAW_LIBS)
	cp $(BUILD_DIR)/libcushaw.so ../../../target/classes/
	rm $(BUILD_DIR)/*.o

clean:
	$(RMR) $(BUILD_DIR)
	$(MAKE) clean -C $(BWA_DIR)
	$(MAKE) clean -C $(CUSHAW_DIR)