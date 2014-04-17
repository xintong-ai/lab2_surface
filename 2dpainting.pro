QT          += opengl


HEADERS     = glwidget.h \
              helper.h \
              window.h \
              mesh.h
SOURCES     = glwidget.cpp \
              helper.cpp \
              main.cpp \
              window.cpp \
              mesh.cpp
			  
win32:INCLUDEPATH += "C:/Program Files/CGAL/include"
win32:INCLUDEPATH += "D:/lib/boost_1_53_0"
win32:INCLUDEPATH += "D:/lib/CGAL-4.2/auxiliary/gmp/include"

win32 {
	

	LIBS += "D:/lib/CGAL-4.2/auxiliary/gmp/lib/libgmp-10.lib"


	debug {
        LIBS += "C:/Program Files/CGAL/lib/CGAL_Core-vc100-mt-gd-4.2.lib"		  
		LIBS += "C:/Program Files/CGAL/lib/CGAL-vc100-mt-gd-4.2.lib"
		LIBS += "D:/lib/boost_1_53_0/stage/lib/libboost_thread-vc100-mt-gd-1_53.lib"
		LIBS += "D:/lib/boost_1_53_0/stage/lib/libboost_system-vc100-mt-gd-1_53.lib"
    }
    release {
		LIBS += "C:/Program Files/CGAL/lib/CGAL_Core-vc100-mt-4.2.lib"		  
		LIBS += "C:/Program Files/CGAL/lib/CGAL-vc100-mt-4.2.lib"
		LIBS += "D:/lib/boost_1_53_0/stage/lib/libboost_thread-vc100-mt-1_53.lib"
		LIBS += "D:/lib/boost_1_53_0/stage/lib/libboost_system-vc100-mt-1_53.lib"
    }
}


			  
# install
#target.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
#sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS 2dpainting.pro
#sources.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
#INSTALLS += target sources

#symbian: include($$QT_SOURCE_TREE/examples/symbianpkgrules.pri)
#maemo5: include($$QT_SOURCE_TREE/examples/maemo5pkgrules.pri)

#symbian: warning(This example might not fully work on Symbian platform)
#simulator: warning(This example might not fully work on Simulator platform)
