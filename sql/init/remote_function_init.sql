CREATE OR REPLACE FUNCTION
  wizepair2.cloudrun.reactor(s string)
  RETURNS string REMOTE
WITH CONNECTION `wizepair2.us.my-connection` OPTIONS ( endpoint = 'https://wizepair2-v1-324998535317.us-east1.run.app',
    max_batching_rows = 1);