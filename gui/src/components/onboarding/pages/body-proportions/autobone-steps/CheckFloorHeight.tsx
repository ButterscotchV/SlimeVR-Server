import {
  ChangeSettingsRequestT,
  HeightRequestT,
  HeightResponseT,
  ModelSettingsT,
  RpcMessage,
  SkeletonHeightT,
} from 'solarxr-protocol';
import { useWebsocketAPI } from '@/hooks/websocket-api';
import { Button } from '@/components/commons/Button';
import { Typography } from '@/components/commons/Typography';
import { Localized, useLocalization } from '@fluent/react';
import { useEffect, useMemo, useState } from 'react';
import { useLocaleConfig } from '@/i18n/config';
import { useHeightContext } from '@/hooks/height';
import { useInterval } from '@/hooks/timeout';
import { TooSmolModal } from './TooSmolModal';
import { MIN_HEIGHT } from '@/components/onboarding/pages/body-proportions/ProportionsChoose';

export function CheckFloorHeightStep({
  nextStep,
  prevStep,
  variant,
}: {
  nextStep: () => void;
  prevStep: () => void;
  variant: 'onboarding' | 'alone';
}) {
  const { l10n } = useLocalization();
  const { floorHeight, hmdHeight, setFloorHeight } = useHeightContext();
  const [fetchedHeight, setFetchedHeight] = useState(false);
  const [fetchHeight, setFetchHeight] = useState(false);
  const { sendRPCPacket, useRPCPacket } = useWebsocketAPI();
  const [isOpen, setOpen] = useState(false);
  const { currentLocales } = useLocaleConfig();

  useEffect(() => setFloorHeight(0), []);

  useInterval(() => {
    if (fetchHeight) {
      sendRPCPacket(RpcMessage.HeightRequest, new HeightRequestT());
    }
  }, 100);

  const mFormat = useMemo(
    () =>
      new Intl.NumberFormat(currentLocales, {
        style: 'unit',
        unit: 'meter',
        maximumFractionDigits: 2,
      }),
    [currentLocales]
  );

  useRPCPacket(RpcMessage.HeightResponse, ({ minHeight }: HeightResponseT) => {
    if (fetchHeight) {
      setFloorHeight((val) =>
        val === null ? minHeight : Math.min(minHeight, val)
      );
    }
  });
  return (
    <>
      <div className="flex flex-col">
        <div className="flex flex-grow">
          <div className="flex flex-grow flex-col gap-4">
            <Typography variant="main-title" bold>
              {l10n.getString(
                'onboarding-automatic_proportions-check_floor_height-title'
              )}
            </Typography>
            <div>
              <Typography color="secondary">
                {l10n.getString(
                  'onboarding-automatic_proportions-check_floor_height-description'
                )}
              </Typography>
              <Localized
                id="onboarding-automatic_proportions-check_floor_height-calculation_warning"
                elems={{ u: <span className="underline"></span> }}
              >
                <Typography color="secondary" bold>
                  Press the button to get your height!
                </Typography>
              </Localized>
            </div>
            <div className="flex flex-grow items-center justify-center">
              <div className="flex flex-col gap-3 items-center">
                {!fetchHeight && (
                  <Button
                    variant="primary"
                    onClick={() => {
                      setFloorHeight(null);
                      setFetchedHeight(false);
                      setFetchHeight(true);
                    }}
                  >
                    <Typography textAlign="text-center">
                      {l10n.getString(
                        fetchedHeight
                          ? 'onboarding-automatic_proportions-check_floor_height-measure-reset'
                          : 'onboarding-automatic_proportions-check_floor_height-measure-start'
                      )}
                    </Typography>
                  </Button>
                )}
                {fetchHeight && (
                  <Button
                    variant="primary"
                    onClick={() => {
                      setFetchHeight(false);
                      setFetchedHeight(true);
                    }}
                  >
                    <Typography textAlign="text-center">
                      {l10n.getString(
                        'onboarding-automatic_proportions-check_floor_height-measure-stop'
                      )}
                    </Typography>
                  </Button>
                )}
                <Typography>
                  {l10n.getString(
                    'onboarding-automatic_proportions-check_floor_height-floor_height'
                  )}
                </Typography>
                <Typography
                  color={fetchHeight ? 'text-status-success' : undefined}
                >
                  {floorHeight === null
                    ? l10n.getString(
                        'onboarding-automatic_proportions-check_height-unknown'
                      )
                    : mFormat.format(floorHeight)}
                </Typography>
              </div>
            </div>
          </div>
          {/* TODO: Get image of person putting controller in floor */}
          {/* <div className="self-center">
            <img
              src="/images/front-standing-pose.webp"
              width={isMobile ? 400 : 300}
              alt="Reset position"
            />
          </div> */}
        </div>

        <div className="flex gap-3 mobile:justify-between">
          <Button
            variant={variant === 'onboarding' ? 'secondary' : 'tertiary'}
            onClick={prevStep}
          >
            {l10n.getString('onboarding-automatic_proportions-prev_step')}
          </Button>
          <Button
            variant={variant === 'onboarding' ? 'secondary' : 'tertiary'}
            disabled={fetchHeight}
            onClick={() => {
              if (!hmdHeight || hmdHeight < MIN_HEIGHT) {
                setOpen(true);
                return;
              }
              setFloorHeight(0);
              const settingsRequest = new ChangeSettingsRequestT();
              settingsRequest.modelSettings = new ModelSettingsT(
                null,
                null,
                null,
                new SkeletonHeightT(hmdHeight, 0)
              );
              sendRPCPacket(RpcMessage.ChangeSettingsRequest, settingsRequest);

              nextStep();
            }}
          >
            {l10n.getString(
              'onboarding-automatic_proportions-check_floor_height-skip_step'
            )}
          </Button>
          <Button
            variant="primary"
            onClick={() => {
              if (
                !hmdHeight ||
                !floorHeight ||
                hmdHeight - floorHeight < MIN_HEIGHT
              ) {
                setOpen(true);
                return;
              }
              const settingsRequest = new ChangeSettingsRequestT();
              settingsRequest.modelSettings = new ModelSettingsT(
                null,
                null,
                null,
                new SkeletonHeightT(hmdHeight, floorHeight)
              );
              sendRPCPacket(RpcMessage.ChangeSettingsRequest, settingsRequest);

              nextStep();
            }}
            disabled={!fetchedHeight && floorHeight == null}
          >
            {l10n.getString(
              'onboarding-automatic_proportions-check_floor_height-next_step'
            )}
          </Button>
        </div>
      </div>
      <TooSmolModal isOpen={isOpen} onClose={() => setOpen(false)} />
    </>
  );
}