if( WITH_ROOT )

    cmessage( STATUS "Creating readdata executable..." )
    add_executable(readData ${CMAKE_CURRENT_SOURCE_DIR}/readData.cpp)
    target_link_libraries(readData ${ROOT_LIBRARIES})
    install( TARGETS readData DESTINATION bin )

    cmessage( STATUS "Creating PerovAna executable..." )
    add_executable(PerovAna ${CMAKE_CURRENT_SOURCE_DIR}/PerovAna.cpp ${CMAKE_CURRENT_SOURCE_DIR}/event.cpp)
    target_link_libraries(PerovAna ${ROOT_LIBRARIES})
    target_link_libraries(PerovAna Spectrum)
    install( TARGETS PerovAna DESTINATION bin )

endif()