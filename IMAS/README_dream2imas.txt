1. Put all the sims to IDS (if necessary), to do that

	python -u ~/DREAM/tools/h5_to_IMAS/export_multiple_sims_to_IMAS.py . --include-folders-file ./simulations_folders2imas.txt --exclude-folders-file exclude_sims.txt --pulse-start 120000 > output_dream2ids 

2. Create manifest files and ingest simulations locally 

	python -u ~/tasks/2026/dream_IMAS/simdb_ingest_multiple_sims.py --uri-map dream_simulation_uris.tsv --template manifest_template.yaml > output_ingest
	You may have to do the followign to delete all simulations
		simdb database clear


3. Validate and check a few sims width simdb commands

4. Push local sim db to remote, CAREFUL, the following pushes all sims!!

	simdb simulation list -l 0 --uuid | awk 'NR > 2 {print $2}' | xargs -r -I{} simdb simulation push {} --username artolaj
