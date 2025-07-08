if( WITH_ROOT )

    cmessage( STATUS "Creating readdata executable..." )
    add_executable(readData ${CMAKE_CURRENT_SOURCE_DIR}/readData.cpp)
    target_link_libraries(readData ${ROOT_LIBRARIES})
    install( TARGETS readData DESTINATION bin )

    cmessage( STATUS "Creating readdata_nch executable..." )
    add_executable(readData_nch ${CMAKE_CURRENT_SOURCE_DIR}/readdata_nch.cpp)
    target_link_libraries(readData_nch ${ROOT_LIBRARIES})
    install( TARGETS readData_nch DESTINATION bin )

    cmessage( STATUS "Creating PerovAna executable..." )
    add_executable(PerovAna ${CMAKE_CURRENT_SOURCE_DIR}/PerovAna.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Event.cpp)
    target_link_libraries(PerovAna ${ROOT_LIBRARIES})
    target_link_libraries(PerovAna Spectrum)
    install( TARGETS PerovAna DESTINATION bin )

    cmessage( STATUS "Creating fft executable..." )
    add_executable(fft ${CMAKE_CURRENT_SOURCE_DIR}/fft.cpp)
    target_link_libraries(fft ${ROOT_LIBRARIES})
    install( TARGETS fft DESTINATION bin )
endif()