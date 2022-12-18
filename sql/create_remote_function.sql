CREATE or replace FUNCTION wizepair2.cloudrun.reactor(s string) RETURNS string
REMOTE WITH CONNECTION `wizepair2.us.my-connection`
OPTIONS (
  endpoint = 'https://wizepair2-v1-46t33xvtyq-ue.a.run.app',
  max_batching_rows = 1);