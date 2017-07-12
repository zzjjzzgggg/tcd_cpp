# estimate
add_library(em_is em_alg.cpp)
target_link_libraries(em_is ${LIBRARIES})
target_compile_definitions(em_is PUBLIC METHOD_IS=1)

add_library(em_us em_alg.cpp)
target_link_libraries(em_us ${LIBRARIES})
target_compile_definitions(em_us PUBLIC METHOD_US=1)

add_library(em_fs em_alg.cpp)
target_link_libraries(em_fs ${LIBRARIES})
target_compile_definitions(em_fs PUBLIC METHOD_FS=1)

add_library(sampler graph_sampler.cpp)
target_link_libraries(sampler ${LIBRARIES})

add_executable(groundtruth main_groundtruth.cpp)
target_link_libraries(groundtruth sampler)

add_executable(estimate_is main_estimate.cpp)
target_link_libraries(estimate_is sampler em_is)

add_executable(estimate_us main_estimate.cpp)
target_link_libraries(estimate_us sampler em_us)

add_executable(estimate_fs main_estimate.cpp)
target_link_libraries(estimate_fs sampler em_fs)

# crlb
add_library(crlb crlb.cpp)
target_link_libraries(crlb ${LIBRARIES})

add_executable(calcrlb_is main_crlb.cpp)
target_link_libraries(calcrlb_is crlb)
target_compile_definitions(calcrlb_is PRIVATE METHOD_IS=1)

add_executable(calcrlb_us main_crlb.cpp)
target_link_libraries(calcrlb_us crlb)
target_compile_definitions(calcrlb_us PRIVATE METHOD_US=1)

add_executable(calcrlb_fs main_crlb.cpp)
target_link_libraries(calcrlb_fs crlb)
target_compile_definitions(calcrlb_fs PRIVATE METHOD_FS=1)

# truncate graph
add_executable(truncate_graph main_truncate_graph.cpp)
target_link_libraries(truncate_graph ${LIBRARIES})

# test
add_executable(test-me main_test.cpp)
target_link_libraries(test-me ${LIBRARIES})
