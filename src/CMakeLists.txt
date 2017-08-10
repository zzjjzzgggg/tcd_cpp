# truncate graph
add_executable(truncate_graph main_truncate_graph.cpp)
target_link_libraries(truncate_graph ${LIBRARIES})

# groundtruth
add_executable(groundtruth main_groundtruth.cpp)
target_link_libraries(groundtruth ${LIBRARIES})

# estimate
add_library(em em_alg.cpp)
target_link_libraries(em ${LIBRARIES})

add_library(em_un em_alg.cpp)
target_link_libraries(em_un ${LIBRARIES})
target_compile_definitions(em_un PUBLIC N_UN=1)

##################################################################################
# n known
## ITS
add_executable(estimate_its main_estimate.cpp)
target_link_libraries(estimate_its em)
target_compile_definitions(estimate_its PRIVATE S_ITS=1)

## ITS-color
add_executable(estimate_itsc main_estimate.cpp)
target_link_libraries(estimate_itsc em)
target_compile_definitions(estimate_itsc PRIVATE S_ITSC=1)

## SGS
add_executable(estimate_sgs main_estimate.cpp)
target_link_libraries(estimate_sgs em)
target_compile_definitions(estimate_sgs PRIVATE S_SGS=1)

##################################################################################
# n unknown
## ITS
add_executable(estimate_un_its main_estimate.cpp)
target_link_libraries(estimate_un_its em_un)
target_compile_definitions(estimate_un_its PRIVATE S_ITS=1)

## ITS-color
add_executable(estimate_un_itsc main_estimate.cpp)
target_link_libraries(estimate_un_itsc em_un)
target_compile_definitions(estimate_un_itsc PRIVATE S_ITSC=1)

## SGS
add_executable(estimate_un_sgs main_estimate.cpp)
target_link_libraries(estimate_un_sgs em_un)
target_compile_definitions(estimate_un_sgs PRIVATE S_SGS=1)

##################################################################################
# crlb
add_executable(calcrlb_its main_crlb.cpp)
target_link_libraries(calcrlb_its ${LIBRARIES})
target_compile_definitions(calcrlb_its PRIVATE S_ITS=1)

add_executable(calcrlb_sgs main_crlb.cpp)
target_link_libraries(calcrlb_sgs ${LIBRARIES})
target_compile_definitions(calcrlb_sgs PRIVATE S_SGS=1)

# test
add_executable(test-me main_test.cpp)
target_link_libraries(test-me ${LIBRARIES})
