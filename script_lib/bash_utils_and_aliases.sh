


func_timing_start () {
	echo "${FUNCNAME[0]} started at $(date)" >> $store/timing
}

func_timing_end () {
	echo "${FUNCNAME[0]} ended at $(date)" >> $store/timing
}
