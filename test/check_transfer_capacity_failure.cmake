# Check both the guard's exit status and its diagnostic. WILL_FAIL accepts
# sanitizer exits, while PASS_REGULAR_EXPRESSION alone ignores the exit code.
if(NOT DEFINED PROGRAM OR NOT DEFINED ARGUMENT OR NOT DEFINED EXPECTED)
    message(FATAL_ERROR "PROGRAM, ARGUMENT and EXPECTED are required")
endif()

execute_process(
    COMMAND "${PROGRAM}" "${ARGUMENT}"
    RESULT_VARIABLE status
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
    TIMEOUT 10)
if(NOT "${status}" STREQUAL "1")
    message(FATAL_ERROR "Expected guard exit 1, got ${status}\n${output}${error}")
endif()
string(REPLACE "\r\n" "\n" error "${error}")
if(NOT "${error}" STREQUAL "${EXPECTED}\n")
    message(FATAL_ERROR "Expected guard diagnostic: ${EXPECTED}\nGot:\n${output}${error}")
endif()
